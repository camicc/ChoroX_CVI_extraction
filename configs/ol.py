"""OphtalmoLaus settings for the original OCT_vol_seg function.

Provide img_path_OL and central_slice_idx_OL separately when calling the
function; neither value was specified in the original analysis snippet.

With dataset="OL", OCT_vol_seg also selects crop_coor_y=(280, 480) and
padding_svg=50 internally. These are not OCT_vol_seg keyword arguments.
"""

OL_CONFIG = {
    "dataset": "OL",
    "nbr_patients": 49,
    "nbr_slices": 5,
    # Pachyvessel segmentation: retain original defaults except C=6.
    "kernel_size_o_pachy": 17,
    "kernel_size_c_pachy": 5,
    "block_size_pachy": 91,
    "cst_sub_mean_pachy": 6,
    # Medium-vessel segmentation: retain original defaults except C=6.
    "kernel_size_o": 5,
    "block_size": 41,
    "cst_sub_mean": 6,
    "kernel_size_o_small": 5,
    "pachy_only": False,
    "medium_only": False,
}
