import cv2
import numpy as np
from pathlib import Path
from typing import Union
from scipy import ndimage
from typing import BinaryIO
from PIL import Image
from lxml import etree
from svglib.svglib import load_svg_file, SvgRenderer
import io
from enum import Enum
from reportlab.graphics import renderPM
from reportlab.graphics import shapes
from contextlib import suppress as suppress_output


def select_central_oct_slices(
    cube_path: Path,
    seg_path: Path,
    nb_slices: int = 5,
):
    """Select only central B-scans from OCT cube folders"""

    if not cube_path.exists() or not seg_path.exists():
        return [], []

    # Count available OCT slices (cheap)
    jpgs = sorted(cube_path.glob("*.jpg"))
    n_slices = len(jpgs)

    if n_slices == 0:
        return [], []

    nb_slices = min(nb_slices, n_slices)

    # Compute central indices
    central_start_idx = (n_slices - nb_slices) // 2
    indices = range(central_start_idx, central_start_idx + nb_slices)

    OCT_slices = []
    seg_slices = []

    for i in indices:
        stem = str(i).zfill(4)
        jpg = cube_path / f"{stem}.jpg"
        svg = seg_path / f"{stem}.svg"

        if not jpg.exists() or not svg.exists():
            continue  # be tolerant

        OCT_slices.append(read_img(jpg))
        seg_slices.append(read_img(svg))

    return OCT_slices, seg_slices


def detect_bm_from_seg(segmentation, smooth_sigma=4):
    """
    Estimate Bruchs Membrane curve from a choroid mask (1=choroid, 0=background).
    """

    # Extract choroidal layer mask
    mask = seg2mask(segmentation, layers='CC_CS') 
    H, W = mask.shape
    big = H * 4
    z = np.arange(H)[:, None]
    # Finds top pixels of the choroidal layer
    top = np.where(mask > 0, z, big).min(axis=0)
    top[top >= big] = np.median(top[top < big])
    # Smooth boundary
    top = ndimage.gaussian_filter1d(top.astype(np.float32), smooth_sigma)
    return top, mask

def flatten_slice(bscan, z_bm, target_row=None, order=1):
    """
    Warp a 2D OCT slice so BM lies flat.
    """
    H, W = bscan.shape
    if target_row is None:
        # Set target BM to median of old coordinates
        target_row = int(np.median(z_bm))

    # Calculate displacement
    yy, xx = np.indices((H, W), dtype=np.float32)
    shift = (z_bm - target_row).astype(np.float32)
    src_y = yy + shift[None, :]
    src_x = xx

    # Shif BM to target row
    flat_OCT = ndimage.map_coordinates(bscan, [src_y, src_x], order=order, mode='reflect')
   
    return flat_OCT

def crop_to_mask(OCT_slice, mask_choro, padding=10, crop_coor_y=None):
    """
    Mask a 2D OCT slice to the choroidal layer and crop it to the mask’s vertical span.

    If (y_min, y_max) crop coordinates are provided, they’re used directly for
    consistency across a volume. Otherwise they’re computed from the mask (with padding).
    """
    # zero out non-choroid pixels
    masked_OCT = OCT_slice * mask_choro

    if crop_coor_y is None:
        # get foreground row indices from mask
        ys, xs = np.where(mask_choro)

        # derive y-bounds (padded)
        y_min, y_max = ys.min() - padding, ys.max() + padding
    else:
        # reuse provided y-bounds (keeps cube consistent)
        y_min, y_max = crop_coor_y

    # crop rows, keep all cols
    masked_OCT = masked_OCT[y_min:y_max, :]

    return masked_OCT


def basic_denoising(volume, k_size_med=5, sigma=1.0): 
    """
    Remove salt & pepper noise followed by slight gaussian smoothing to have more well-defined edges.
    """
    # median filter removes salt and pepper noise
    m = ndimage.median_filter(volume, size=k_size_med)

    # Slight gaussian filter (anisotropic so less weight in the z direction)
    g = ndimage.gaussian_filter(m, sigma=(sigma/2, sigma, sigma))

    return g


def depth_equalize_cube(volume_lin, mask_valid=None,
                        pct=60, smooth_sigma=8, gain_clip=(0.5, 3.0)):
    """
    Intensity depth-equalize a 3-D OCT cube (N, H, W) along the depth axis for a region of interest (choroidal layer).
    """
    v = np.asarray(volume_lin, np.float32)
    N, H, W = v.shape

    if mask_valid is None:
        mask_valid = np.ones_like(v, bool)
    else:
        mask_valid = np.array(mask_valid).astype(bool)

    # Compute per-depth intensity envelope
    E = np.zeros(H, np.float32)
    for z in range(H):
        # gather all valid voxels at this depth
        vals = v[:, z, :][mask_valid[:, z, :]]
        if vals.size < 64:
            E[z] = np.median(v[:, z, :])
        else:
            E[z] = np.percentile(vals, pct)

    # Smooth envelope along depth
    E = ndimage.gaussian_filter1d(E, smooth_sigma)

    # Compute gain to flatten depth response
    target = np.median(E[E > 0])
    g = np.clip(target / (E + 1e-8), *gain_clip).astype(np.float32)

    # Apply gain to cube
    equalized = np.clip(v * g[None, :, None], 0, None)

    return equalized, E, g

def preprocess_oct(
    oct_cube=None, oct_seg=None,
    bm_layer_height=300, crop_coor_y=(280, 380),
    remove_salt_pepper=True, depth_contrast_enhance=True, adjust_contrast=True,
    smooth_sigma_bm=4, order_flatten=1, padding_crop=10, k_size_med=5, sigma_denoise=1.0,
    pct_equalize=60, smooth_sigma_equalize=8, gain_clip_equalize=(0.5, 3.0)):

    flat_bscans = []
    flat_masks = []
    # Align B-scans w.r.t BM & Crop to Choroidal Layer
    for bscan, seg in zip(oct_cube, oct_seg):
        # Detect BM using upper boundary of choroid layer
        bm_y, mask_choro = detect_bm_from_seg(seg, smooth_sigma=smooth_sigma_bm)
        # Flatten b-scan & mask w.r.t its BM
        flat_bscan = flatten_slice(bscan, bm_y, target_row=bm_layer_height, order=order_flatten)
        flat_mask  = flatten_slice(mask_choro, bm_y, target_row=bm_layer_height, order=0)
        # Keep choroidal layer only and crop to mask 
        flat_bscan = crop_to_mask(flat_bscan, flat_mask, padding=padding_crop, crop_coor_y=crop_coor_y)
        flat_mask  = flat_mask[crop_coor_y[0]:crop_coor_y[1], :]

        # Add to lists
        flat_bscans.append(flat_bscan)
        flat_masks.append(flat_mask)

    # Stack to get flattened 3D volume
    flat_volume = np.stack(flat_bscans, axis=0)
    flat_mask_volume = np.stack(flat_masks, axis=0)

    # Denoising
    denoised_vol = flat_volume

    # Remove salt & pepper noise
    if remove_salt_pepper:
        denoised_vol = basic_denoising(denoised_vol, k_size_med=k_size_med, sigma=sigma_denoise*2)

    # Depth equalization
    if depth_contrast_enhance:
        denoised_vol, _, _ = depth_equalize_cube(
            denoised_vol, mask_valid=flat_mask_volume,
            pct=pct_equalize, smooth_sigma=smooth_sigma_equalize,
            gain_clip=gain_clip_equalize
        )

    if adjust_contrast:
        # Stretch cube's intensity range to (0, 255) to increase contrast
        denoised_vol = denoised_vol.astype(np.float32)
        denoised_vol = (denoised_vol - denoised_vol.min()) / (denoised_vol.max() - denoised_vol.min()) * 255
        denoised_vol = denoised_vol.astype(np.uint8)

    # Slight gaussian filter (anisotropic so less weight in the z direction)
    denoised_vol = ndimage.gaussian_filter(denoised_vol, sigma=(1, 2, 2))

    return denoised_vol, flat_mask_volume

def segment_oct(oct_cube=None, oct_seg=None, 
                # Pachy vessels detection parameters
                kernel_size_o_pachy=17, kernel_size_c_pachy=5, block_size_pachy=91, cst_sub_mean_pachy=13,

                # Medium vessels detection parameters
                kernel_size_o=5, block_size=41, cst_sub_mean=16, kernel_size_o_small=5,

                # Modes
                pachy_only=False, medium_only=False):

    """ Perform preprocessing, thresholding and morphological operators on processed OCT cubes"""
    # Preprocess OCT cubes
    processed_cube, flat_mask_volume = preprocess_oct(oct_cube, oct_seg)

    thresholded_cube = []

    for i, slice in enumerate(processed_cube):
        # Ensure grayscale and 8-bit
        img_8bit = cv2.normalize(slice, None, 0, 255, cv2.NORM_MINMAX)  # scale to 0–255
        img_8bit = img_8bit.astype(np.uint8)
        img_8bit = cv2.medianBlur(img_8bit, 3)

        
        # Threshold slice
        thresh_pachy = cv2.adaptiveThreshold(
            img_8bit,
            255,                           
            cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
            cv2.THRESH_BINARY_INV,         
            blockSize=block_size_pachy,                
            C=cst_sub_mean_pachy                     
        )

        # Threshold slice
        thresh = cv2.adaptiveThreshold(
            img_8bit,
            255,                           
            cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
            cv2.THRESH_BINARY_INV,         
            blockSize=block_size,                
            C=cst_sub_mean                     
        )

        
        # Create kernels with ellipses to ressemble circular vessel shapes
        kernel_closing_pachy = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size_c_pachy, kernel_size_c_pachy))
        kernel_opening_pachy = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size_o_pachy, kernel_size_o_pachy))
        
        kernel_opening = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size_o, kernel_size_o))

        kernel_opening_small = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size_o_small, kernel_size_o_small))

        # Pachy vessels: Unlink connected circles and remove medium elements 
        thresh_pachy = cv2.bitwise_and(thresh_pachy, (flat_mask_volume[i]).astype(np.uint8))
        closed  = cv2.morphologyEx(thresh_pachy , cv2.MORPH_CLOSE, kernel_closing_pachy)
        pachy  = cv2.morphologyEx(closed, cv2.MORPH_OPEN, kernel_opening_pachy)

        # Medium vessels: Unlink connected circles and remove medium element
        medium  = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel_opening)
        medium = cv2.bitwise_and(medium, cv2.bitwise_not(pachy))
        medium = cv2.morphologyEx(medium, cv2.MORPH_OPEN, kernel_opening_small)

        # Compute pachy vessels only
        if pachy_only:
            thresholded_cube.append(pachy)
            continue

        # Compute small and medium vessels only
        if medium_only:
            thresholded_cube.append(medium)
            continue
            
       
        # Dilate pachy slightly to define regions that should be empty between pachy and medium vessels
        no_touch_zone = cv2.dilate(pachy, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)))

        # Remove medium regions that intersect pachy
        medium = cv2.bitwise_and(medium, cv2.bitwise_not(no_touch_zone))
            
        all_vessels = pachy + medium 

               
        # Remove boundary artefacts
        # Keep only pixels farther than N from the boundar
        dist = cv2.distanceTransform(flat_mask_volume[i].astype(np.uint8), cv2.DIST_L2, 5)
        enlarged_mask = (dist > 8).astype(np.uint8)
        all_vessels = cv2.bitwise_and(all_vessels, (enlarged_mask).astype(np.uint8))

        # Deconnect vessels
        all_vessels = cv2.morphologyEx(all_vessels, cv2.MORPH_OPEN, kernel_opening_small)

        thresholded_cube.append(all_vessels)

    return thresholded_cube, flat_mask_volume

def calculate_cvi(segmented_cube, cropped_flat_masks):
    """
    Calculate choroidal vascularity index (CVI) per x coordinate.

    For each B-scan, CVI(x) = vessel_pixels_in_column / mask_thickness(x),
    """
    # one array per B-scan
    cvi_mean_bscans = []
    cvi_std_bscans = []

    for segmented_vessels, mask in zip(segmented_cube, cropped_flat_masks):
        seg = np.array(segmented_vessels)
        msk = np.array(mask)

        # Handle empty masks
        if msk.sum() == 0:
            continue

        # Sum over the vertical (y) axis to get per-column counts
        vessel_per_x = np.sum(seg, axis=0)   # number of vessel pixels in each x-column
        mask_thickness = np.sum(msk, axis=0) # total choroid pixels in each x-column

        # Avoid divide-by-zero
        with np.errstate(divide='ignore', invalid='ignore'):
            cvi_per_x = np.where(mask_thickness > 0, vessel_per_x / mask_thickness, 0)
        
        cvi_mean_bscans.append(np.mean(cvi_per_x))
        cvi_std_bscans.append(np.std(cvi_per_x))

    cvi_mean = np.median(cvi_mean_bscans) if cvi_mean_bscans else None
    cvi_std = np.median(cvi_std_bscans) if cvi_std_bscans else None

    return cvi_mean, cvi_std


def calculate_pachy_area(segmented_cube, spacings=None):
    """ Calculate pachy area by finding the largest segmented and summing its pixels"""
    pachy_areas = []
    for segmented_oct in segmented_cube:
        # Find all connected component
        numLabels, labels, stats, centroids = cv2.connectedComponentsWithStats(segmented_oct, connectivity=8, ltype=cv2.CV_16U)

        # Get area of largest component
        areas = stats[1:, cv2.CC_STAT_AREA]  # skip background (label 0)

        if areas.size == 0:
            continue
        else:
            # Convert area to mm^2 (z:height, x:width) 
            # Spacings are in [mm/pixel] 
            pachy_areas.append(areas.max() * spacings[2] * spacings[3])

    # Take the median of over bscans
    pachy_area = np.median(pachy_areas) if len(pachy_areas) else None
    
    return pachy_area

def read_img(file: Union[Path, str, BinaryIO]) -> Union[np.ndarray, etree._Element]:
    """Reads an image from a file or a binary IO object."""

    # Convert str to Path object
    if isinstance(file, str):
        file = Path(file)

    # Read a file
    if isinstance(file, Path):
        if file.suffix.lower() in ['.svg', '.svgz']:
            svg = load_svg_file(str(file))
            for p in svg[1]:
                p.attrib['fill-opacity'] = '1'
            return svg
        else:
            img = Image.open(file)
            img = np.asarray(img)
            return img

    # Read from a file-like object
    elif all([
            hasattr(file, 'read'),
            hasattr(file, 'seek'),
    ]):
        file.seek(0)
        img = Image.open(file)
        img = np.asarray(img)
        return img
    
def seg2mask(svg: etree._Element, layers: Union[str, list[str], None] = None) -> np.ndarray:
    """Returns the mask of the specified layers in a segmentation.

    Args:
        svg: Segmentation svg.
        layers: The layer(s) to be masked.
            If ``None``, it will return all the layers except the background. Defaults to ``None`` .

    Returns:
        The mask array of the specified layers.

    .. seealso::
        `mask2bnd <src.cohortextractor.utils.imageprocessing.mask2bnd>`
            Returns the boundary of a mask.
        `bnd2curve <src.cohortextractor.utils.imageprocessing.mask2bnd>`
            Extracts the indices of the boundary pixels.
    """

    # Interpret the layers
    if not layers:
        layers = [layer for layer in RetinalLayer.get_names() if layer != 'BG']
    elif isinstance(layers, str):
        layers = [layers]
    # Check layer names
    assert all(layer in RetinalLayer.get_names() for layer in layers),\
        f'Layers should be a subset of: {RetinalLayer.get_names()}'
    # Convert inputs
    layer_tags = [RetinalLayer[layer].value['tag'] for layer in layers]
    # Make other layers transparent
    for p in svg[1]:
        if p.attrib['class'] not in layer_tags:
            p.attrib['fill-opacity'] = '0'  # transparent
        if p.attrib['class'] in layer_tags:
            p.attrib['fill-opacity'] = '1'  # opaque
    svgRenderer = SvgRenderer('')
    drawing = svgRenderer.render(svg)
    img = draw2img(drawing, fmt='PNG')
    mask = np.array(Image.open(img))
    mask = np.round(np.amax(1 - mask/255, axis=2))
    mask = np.array(mask, dtype=bool)
    return mask


def draw2img(drawing: shapes.Drawing, **kwargs) -> io.BytesIO:
    """
    Creates a pixmap and draw drawing to it then return it as a BytesIO.

    .. seealso::
        `Drawing <https://pyng.tech/docs/reportlab/graphics/shapes.m.html#reportlab.graphics.shapes.Drawing>`_
            Drawing objects in the reportlab documentation.
        `renderPM.drawToFile <https://pyng.tech/docs/reportlab/graphics/renderPM.m.html#reportlab.graphics.renderPM.drawToFile>`_
            Base method in the reportlab documentation.
    """

    img = io.BytesIO()
    with suppress_output():
        renderPM.drawToFile(drawing, img, **kwargs)

    return img


class RetinalLayer(Enum):
    """
    The title, the SVG tag (in segmentations), and the type of retinal layers.
    """

    # NOTE: The order is important.
    BG = {
        'title': 'Background',
        'tag': 'layer layer-0',
        'type': None,
        }
    RNFL = {
        'title': 'Retinal Nerve Fibre Layer',
        'tag': 'layer layer-1',
        'type': 'layer',
        }
    GCL_IPL = {
        'title': 'Ganglion Cell Layer and Inner Plexiform Layer',
        'tag': 'layer layer-2',
        'type': 'layer',
        }
    INL_OPL = {
        'title': 'Inner Nuclear Layer and Outer Plexiform Layer',
        'tag': 'layer layer-3',
        'type': 'layer',
        }
    ONL = {
        'title': 'Outer Nuclear Layer',
        'tag': 'layer layer-4',
        'type': 'layer',
        }
    PR_RPE = {
        'title': 'Photoreceptors and Retinal Pigment Epithelium',
        'tag': 'layer layer-5',
        'type': 'layer',
        }
    CC_CS = {
        'title': 'Choriocapillaris and Choroidal Stroma',
        'tag': 'layer layer-6',
        'type': 'layer',
        }
    IRF = {
        'title': 'Intraretinal Fluid',
        'tag': 'layer layer-7',
        'type': 'fluid',
        }
    SRF = {
        'title': 'Subretinal Fluid',
        'tag': 'layer layer-8 fluid',
        'type': 'fluid',
        }
    PED = {
        'title': 'Pigment Epithelial Detachment',
        'tag': 'layer layer-9 fluid',
        'type': 'fluid',
        }

    # CHECK: What is this layer?
    # NOTE: Extra number in the thickness or volume json files that
    # does not match with anything on Discovery UI
    Unknown = {
        'title': None,
        'tag': None,
        'type': None,
        }

    @classmethod
    def get_names(cls, plus: bool = False) -> list[str]:
        """
        Returns a list of the names.
        If ``plus`` is True, the underscore in the names will be replaced by ``"+"``.
        """

        return [layer.name.replace('_', '+') if plus else layer.name for layer in cls]

    @classmethod
    def get_names_layers(cls, plus: bool = False) -> list[str]:
        """
        Returns a list of the names (only layers).
        If ``plus`` is True, the underscore in the names will be replaced by ``"+"``.
        """

        return [layer.name.replace('_', '+') if plus else layer.name for layer in cls if layer.value['type'] == 'layer']

    @classmethod
    def get_names_fluids(cls, plus: bool = False) -> list[str]:
        """
        Returns a list of the names (only fluids).
        If ``plus`` is True, the underscore in the names will be replaced by ``"+"``.
        """

        return [layer.name.replace('_', '+') if plus else layer.name for layer in cls if layer.value['type'] == 'fluid']

    @classmethod
    def get_titles(cls) -> list[str]:
        """Returns a list of the titles."""

        return [layer.value['title'] for layer in cls]

    @classmethod
    def get_tags(cls) -> list[str]:
        """Returns a list of the tags."""

        return [layer.value['tag'] for layer in cls]
