"""Original CVI image-processing operations, exposed through a small public API.

The numerical steps intentionally follow the original implementation. In
particular, do not change the bitwise mask operations, normalization, cropping,
interpolation, or the order of morphological filters without a regression test.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from io import BytesIO
from pathlib import Path
from typing import BinaryIO, Sequence

import cv2
import matplotlib.pyplot as plt
import numpy as np
from lxml import etree
from PIL import Image
from reportlab.graphics import renderPM
from scipy import ndimage
from scipy.ndimage import gaussian_filter, gaussian_filter1d, map_coordinates, median_filter, zoom
from svglib.svglib import SvgRenderer, load_svg_file
from reportlab import rl_config
rl_config.renderPMBackend = "_renderPM"


class RetinalLayer(Enum):
    """Discovery SVG layer names and class tags, in the original order."""

    BG = {"title": "Background", "tag": "layer layer-0", "type": None}
    RNFL = {"title": "Retinal Nerve Fibre Layer", "tag": "layer layer-1", "type": "layer"}
    GCL_IPL = {
        "title": "Ganglion Cell Layer and Inner Plexiform Layer",
        "tag": "layer layer-2", "type": "layer",
    }
    INL_OPL = {
        "title": "Inner Nuclear Layer and Outer Plexiform Layer",
        "tag": "layer layer-3", "type": "layer",
    }
    ONL = {"title": "Outer Nuclear Layer", "tag": "layer layer-4", "type": "layer"}
    PR_RPE = {
        "title": "Photoreceptors and Retinal Pigment Epithelium",
        "tag": "layer layer-5", "type": "layer",
    }
    CC_CS = {
        "title": "Choriocapillaris and Choroidal Stroma",
        "tag": "layer layer-6", "type": "layer",
    }
    IRF = {"title": "Intraretinal Fluid", "tag": "layer layer-7", "type": "fluid"}
    SRF = {"title": "Subretinal Fluid", "tag": "layer layer-8 fluid", "type": "fluid"}
    PED = {"title": "Pigment Epithelial Detachment", "tag": "layer layer-9 fluid", "type": "fluid"}
    Unknown = {"title": None, "tag": None, "type": None}

    @classmethod
    def get_names(cls, plus: bool = False) -> list[str]:
        return [layer.name.replace("_", "+") if plus else layer.name for layer in cls]

    @classmethod
    def get_names_layers(cls, plus: bool = False) -> list[str]:
        return [
            layer.name.replace("_", "+") if plus else layer.name
            for layer in cls if layer.value["type"] == "layer"
        ]

    @classmethod
    def get_names_fluids(cls, plus: bool = False) -> list[str]:
        return [
            layer.name.replace("_", "+") if plus else layer.name
            for layer in cls if layer.value["type"] == "fluid"
        ]

    @classmethod
    def get_titles(cls) -> list[str | None]:
        return [layer.value["title"] for layer in cls]

    @classmethod
    def get_tags(cls) -> list[str | None]:
        return [layer.value["tag"] for layer in cls]


def read_img(file: Path | str | BinaryIO) -> np.ndarray | etree._Element:
    """Read an OCT image or its SVG segmentation, as in the original pipeline."""
    if isinstance(file, (str, Path)):
        file = Path(file)
        if file.suffix.lower() in (".svg", ".svgz"):
            svg = load_svg_file(str(file))
            for path in svg[1]:
                path.attrib["fill-opacity"] = "1"
            return svg
        with Image.open(file) as img:
            return np.asarray(img)

    if hasattr(file, "read") and hasattr(file, "seek"):
        file.seek(0)
        with Image.open(file) as img:
            return np.asarray(img)

    raise TypeError("Expected an image path or a seekable binary file")


def seg2mask(svg: etree._Element, layers: str | list[str] | None = None) -> np.ndarray:
    """Rasterize selected SVG layers without modifying the supplied SVG."""
    svg = deepcopy(svg)
    if not layers:
        layers = [name for name in RetinalLayer.get_names() if name != "BG"]
    elif isinstance(layers, str):
        layers = [layers]

    assert all(layer in RetinalLayer.get_names() for layer in layers), (
        f"Layers should be a subset of: {RetinalLayer.get_names()}"
    )
    layer_tags = [RetinalLayer[layer].value["tag"] for layer in layers]

    for path in svg[1]:
        path.attrib["fill-opacity"] = "1" if path.attrib["class"] in layer_tags else "0"

    drawing = SvgRenderer("").render(svg)
    png = BytesIO()
    renderPM.drawToFile(drawing, png, fmt="PNG")
    png.seek(0)
    with Image.open(png) as img:
        mask = np.array(img)

    # Preserve the original RGB raster-to-boolean conversion exactly.
    mask = np.round(np.amax(1 - mask / 255, axis=2))
    return np.array(mask, dtype=bool)


def detect_bm_from_seg(segmentation: etree._Element, smooth_sigma: float = 4):
    """Estimate Bruch's membrane using the upper edge of the choroid mask."""
    mask = seg2mask(segmentation, layers="CC_CS")
    height, _ = mask.shape
    big = height * 4
    z = np.arange(height)[:, None]
    top = np.where(mask > 0, z, big).min(axis=0)
    top[top >= big] = np.median(top[top < big])
    top = gaussian_filter1d(top.astype(np.float32), smooth_sigma)
    return top, mask


def flatten_slice(bscan: np.ndarray, z_bm: np.ndarray, target_row: int | None = None, order: int = 1):
    """Flatten a B-scan relative to Bruch's membrane."""
    height, width = bscan.shape
    if target_row is None:
        target_row = int(np.median(z_bm))
    yy, xx = np.indices((height, width), dtype=np.float32)
    shift = (z_bm - target_row).astype(np.float32)
    src_y = yy + shift[None, :]
    return map_coordinates(bscan, [src_y, xx], order=order, mode="reflect")


def crop_to_mask(oct_slice, mask_choro, padding_y=10, padding_x=25, crop_coor_y=None):
    """Mask to the choroid and crop the slice using the original coordinates."""
    masked_oct = oct_slice * mask_choro
    if crop_coor_y is None:
        ys, _ = np.where(mask_choro)
        y_min, y_max = ys.min() - padding_y, ys.max() + padding_y
    else:
        y_min, y_max = crop_coor_y
    y_min = max(0, y_min)
    y_max = min(mask_choro.shape[0], y_max)
    masked_oct = masked_oct[y_min:y_max, padding_x:-padding_x]
    mask_choro = mask_choro[y_min:y_max, padding_x:-padding_x]
    return masked_oct, mask_choro


def basic_denoising(volume, k_size_med=5, sigma=1.0):
    """Median filtering followed by anisotropic Gaussian filtering."""
    median = median_filter(volume, size=k_size_med)
    return gaussian_filter(median, sigma=(sigma / 2, sigma, sigma))


def depth_equalize_cube(volume_lin, mask_valid=None, pct=60, smooth_sigma=8, gain_clip=(0.5, 3.0)):
    """Apply the original per-depth intensity gain to a 3-D OCT cube."""
    v = np.asarray(volume_lin, np.float32)
    n, height, width = v.shape
    if mask_valid is None:
        mask_valid = np.ones_like(v, bool)
    else:
        mask_valid = np.array(mask_valid).astype(bool)

    envelope = np.zeros(height, np.float32)
    for z in range(height):
        vals = v[:, z, :][mask_valid[:, z, :]]
        if vals.size < 64:
            envelope[z] = np.median(v[:, z, :])
        else:
            envelope[z] = np.percentile(vals, pct)
    envelope = gaussian_filter1d(envelope, smooth_sigma)
    target = np.median(envelope[envelope > 0])
    gain = np.clip(target / (envelope + 1e-8), *gain_clip).astype(np.float32)
    equalized = np.clip(v * gain[None, :, None], 0, None)
    return equalized, envelope, gain


def _preprocess_oct_full(
    oct_cube,
    oct_seg,
    bm_layer_height=300,
    crop_coor_y=(280, 380),
    remove_salt_pepper=True,
    depth_contrast_enhance=True,
    adjust_contrast=True,
    smooth_sigma_bm=4,
    order_flatten=1,
    padding_crop=10,
    k_size_med=5,
    sigma_denoise=1.0,
    pct_equalize=60,
    smooth_sigma_equalize=8,
    gain_clip_equalize=(0.5, 3.0),
    padding_svg=25,
):
    """Full original return: processed cube, choroid masks, pre-denoising scans."""
    flat_bscans = []
    flat_masks = []
    for bscan, seg in zip(oct_cube, oct_seg):
        bm_y, mask_choro = detect_bm_from_seg(seg, smooth_sigma=smooth_sigma_bm)
        flat_bscan = flatten_slice(bscan, bm_y, target_row=bm_layer_height, order=order_flatten)
        flat_mask = flatten_slice(mask_choro, bm_y, target_row=bm_layer_height, order=0)
        flat_bscan, flat_mask = crop_to_mask(
            flat_bscan,
            flat_mask,
            padding_y=padding_crop,
            padding_x=padding_svg,
            crop_coor_y=crop_coor_y,
        )
        flat_bscans.append(flat_bscan)
        flat_masks.append(flat_mask)

    flat_volume = np.stack(flat_bscans, axis=0)
    flat_mask_volume = np.stack(flat_masks, axis=0)
    denoised_vol = flat_volume
    if remove_salt_pepper:
        denoised_vol = basic_denoising(denoised_vol, k_size_med=k_size_med, sigma=sigma_denoise * 2)
    if depth_contrast_enhance:
        denoised_vol, _, _ = depth_equalize_cube(
            denoised_vol,
            mask_valid=flat_mask_volume,
            pct=pct_equalize,
            smooth_sigma=smooth_sigma_equalize,
            gain_clip=gain_clip_equalize,
        )
    if adjust_contrast:
        denoised_vol = denoised_vol.astype(np.float32)
        denoised_vol = (denoised_vol - denoised_vol.min()) / (denoised_vol.max() - denoised_vol.min()) * 255
        denoised_vol = denoised_vol.astype(np.uint8)
    denoised_vol = ndimage.gaussian_filter(denoised_vol, sigma=(1, 2, 2))
    return denoised_vol, flat_mask_volume, flat_bscans


def preprocess_oct(oct_cube, oct_seg, *, crop_coor_y=(280, 380), padding_svg=1, **kwargs):
    """Return processed OCT and masks with the same crop as ``segment_oct``.

    The *original* stand-alone default was padding_svg=25; segmentation used 1.
    This public wrapper intentionally defaults to 1 for matching debug images.
    All numerical processing steps remain unchanged.
    """
    processed, masks, _ = _preprocess_oct_full(
        oct_cube, oct_seg, crop_coor_y=crop_coor_y, padding_svg=padding_svg, **kwargs
    )
    return processed, masks


def segment_oct(
    oct_cube,
    oct_seg,
    kernel_size_o_pachy=17,
    kernel_size_c_pachy=5,
    block_size_pachy=91,
    cst_sub_mean_pachy=13,
    kernel_size_o=5,
    block_size=41,
    cst_sub_mean=16,
    kernel_size_o_small=5,
    pachy_only=False,
    medium_only=False,
    crop_coor_y=(280, 380),
    padding_svg=1,
    *,
    return_intermediates=False,
):
    """Segment choroidal vessels with the original morphology and thresholds.

    Returns (vessel_slices, mask_volume) by default. With
    ``return_intermediates=True``, appends the processed and flattened volumes.
    Values are preserved as generated by the original OpenCV bitwise operations.
    """
    processed_cube, flat_mask_volume, flat_volume = _preprocess_oct_full(
        oct_cube, oct_seg, crop_coor_y=crop_coor_y, padding_svg=padding_svg
    )
    thresholded_cube = []
    for i, oct_slice in enumerate(processed_cube):
        img_8bit = cv2.normalize(oct_slice, None, 0, 255, cv2.NORM_MINMAX)
        img_8bit = img_8bit.astype(np.uint8)
        img_8bit = cv2.medianBlur(img_8bit, 3)

        thresh_pachy = cv2.adaptiveThreshold(
            img_8bit, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            cv2.THRESH_BINARY_INV, blockSize=block_size_pachy, C=cst_sub_mean_pachy,
        )
        thresh = cv2.adaptiveThreshold(
            img_8bit, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            cv2.THRESH_BINARY_INV, blockSize=block_size, C=cst_sub_mean,
        )

        closing_pachy = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (kernel_size_c_pachy, kernel_size_c_pachy)
        )
        opening_pachy = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (kernel_size_o_pachy, kernel_size_o_pachy)
        )
        opening = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size_o, kernel_size_o))
        opening_small = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (kernel_size_o_small, kernel_size_o_small)
        )

        # Keep the *bitwise* uint8 operations: replacing these with boolean
        # masks changes numerical results (0/1 versus 0/255).
        thresh_pachy = cv2.bitwise_and(thresh_pachy, flat_mask_volume[i].astype(np.uint8))
        closed = cv2.morphologyEx(thresh_pachy, cv2.MORPH_CLOSE, closing_pachy)
        pachy = cv2.morphologyEx(closed, cv2.MORPH_OPEN, opening_pachy)

        medium = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, opening)
        medium = cv2.bitwise_and(medium, cv2.bitwise_not(pachy))
        medium = cv2.morphologyEx(medium, cv2.MORPH_OPEN, opening_small)
        if pachy_only:
            thresholded_cube.append(pachy)
            continue
        if medium_only:
            thresholded_cube.append(medium)
            continue

        no_touch_zone = cv2.dilate(
            pachy, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
        )
        medium = cv2.bitwise_and(medium, cv2.bitwise_not(no_touch_zone))
        all_vessels = pachy + medium

        dist = cv2.distanceTransform(flat_mask_volume[i].astype(np.uint8), cv2.DIST_L2, 5)
        enlarged_mask = (dist > 8).astype(np.uint8)
        all_vessels = cv2.bitwise_and(all_vessels, enlarged_mask.astype(np.uint8))
        all_vessels = cv2.morphologyEx(all_vessels, cv2.MORPH_OPEN, opening_small)
        thresholded_cube.append(all_vessels)

    if return_intermediates:
        return thresholded_cube, flat_mask_volume, processed_cube, flat_volume
    return thresholded_cube, flat_mask_volume


def calculate_cvi(segmented_cube, cropped_flat_masks):
    """Median of per-slice CVI(x) means and within-slice CVI(x) SDs."""
    cvi_mean_bscans = []
    cvi_std_bscans = []
    for segmented_vessels, mask in zip(segmented_cube, cropped_flat_masks):
        seg = np.array(segmented_vessels)
        msk = np.array(mask)
        if msk.sum() == 0:
            cvi_mean_bscans.append(None)
            cvi_std_bscans.append(None)
            continue
        vessel_per_x = np.sum(seg, axis=0)
        mask_thickness = np.sum(msk, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            cvi_per_x = np.where(mask_thickness > 0, vessel_per_x / mask_thickness, 0)
        cvi_mean_bscans.append(np.mean(cvi_per_x))
        cvi_std_bscans.append(np.std(cvi_per_x))

    valid_mean = [value for value in cvi_mean_bscans if value is not None]
    valid_std = [value for value in cvi_std_bscans if value is not None]
    return (
        np.median(valid_mean) if valid_mean else None,
        np.median(valid_std) if valid_std else None,
    )


def calculate_pachy_area(segmented_cube, spacings: Sequence[float]):
    """Legacy connected-component areas (µm²), including background.

    WARNING: This deliberately retains the original ``stats[:, AREA].max()``.
    For a true largest *vessel* area, background label 0 must be excluded; that
    scientific correction would change the original results.
    """
    pachy_areas = []
    for segmented_oct in segmented_cube:
        _, _, stats, _ = cv2.connectedComponentsWithStats(
            segmented_oct, connectivity=8, ltype=cv2.CV_16U
        )
        areas = stats[:, cv2.CC_STAT_AREA]
        pachy_areas.append(areas.max() * spacings[2] * spacings[3] * 1e6)
    return pachy_areas


@dataclass(frozen=True)
class CVIResult:
    """Summary of one OCT volume (legacy area includes background)."""

    cvi_mean: float | None
    cvi_spatial_sd: float | None
    pachy_area_mm2: float
    pachy_areas_mm2: tuple[float, ...]


def extract_cvi(bscans, segmentations, *, spacings: Sequence[float]) -> CVIResult:
    """Segment B-scans and return the original CVI metrics.

    The legacy area function returns one area per slice in µm². Its median is
    converted to mm² for the scalar field; the original did not aggregate it.
    """
    vessels, masks = segment_oct(bscans, segmentations)
    cvi_mean, cvi_spatial_sd = calculate_cvi(vessels, masks)
    areas_mm2 = np.asarray(calculate_pachy_area(vessels, spacings), dtype=float) / 1e6
    return CVIResult(
        cvi_mean=cvi_mean,
        cvi_spatial_sd=cvi_spatial_sd,
        pachy_area_mm2=float(np.median(areas_mm2)),
        pachy_areas_mm2=tuple(float(area) for area in areas_mm2),
    )


def save_segmentation_overlays(processed_cube, segmented_cube, flat_mask_volume, output_dir, *, scan_indices=None):
    """Save diagnostic vessel masks and OCT overlays; no impact on biomarkers."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if scan_indices is None:
        scan_indices = range(len(segmented_cube))
    for i, oct_slice, vessels, mask in zip(
        scan_indices, processed_cube, segmented_cube, flat_mask_volume
    ):
        grayscale = cv2.normalize(oct_slice, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        rgb = cv2.cvtColor(grayscale, cv2.COLOR_GRAY2RGB)
        inside = np.asarray(mask, dtype=bool)
        vessel_pixels = (np.asarray(vessels) > 0) & inside

        # Choroid light blue, vessel pixels red, outside original grayscale.
        overlay = rgb.copy()
        overlay[inside] = np.round(0.8 * rgb[inside] + 0.2 * np.array([0, 190, 230])).astype(np.uint8)
        overlay[vessel_pixels] = [255, 45, 45]
        Image.fromarray(overlay).save(output_dir / f"{i:04d}_overlay.png")
        Image.fromarray((np.asarray(vessels) > 0).astype(np.uint8) * 255).save(
            output_dir / f"{i:04d}_vessels.png"
        )


def plot_enface_oct(segmented_cube, patient_id, save_path=None, spacings=(0.046875, 0, 0.0035, 0.012), *, show=False):
    """Mean-depth en face projection, with the original linear z-resampling."""
    oct_volume = np.array(segmented_cube)
    _, _, dz, dx = spacings
    zoom_factors = [dx / dz, 1.0, 1.0]
    iso_volume = zoom(oct_volume, zoom=zoom_factors, order=1)
    enface_img = np.mean(iso_volume, axis=1)
    if save_path is not None:
        plt.imsave(save_path, enface_img, cmap="gray")
    if show:
        plt.figure(figsize=(7, 7))
        plt.imshow(enface_img, cmap="gray")
        plt.title(f"Enface projection: {patient_id}")
        plt.axis("off")
        plt.show()
    return enface_img