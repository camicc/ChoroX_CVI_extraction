"""Generate five artificial OCT B-scans and matching choroid SVG masks.

These are synthetic illustrations, not OCT measurements or clinical images.
Run: python examples/make_synthetic_example.py
"""

from pathlib import Path
import json

import numpy as np
from PIL import Image


EXAMPLE_DIR = Path(__file__).resolve().parent / "synthetic_case"
SVG = '''<svg xmlns="http://www.w3.org/2000/svg" width="512" height="512" viewBox="0 0 512 512">
  <rect width="512" height="512" fill="white" />
  <g>
    <path class="layer layer-6" d="M 0 300 L 512 300 L 512 390 L 0 390 Z" fill="black" />
  </g>
</svg>
'''


def generate_example(directory: Path = EXAMPLE_DIR) -> None:
    cube = directory / "oct_volume"
    seg = directory / "svg_segmentation"
    cube.mkdir(parents=True, exist_ok=True)
    seg.mkdir(parents=True, exist_ok=True)
    yy, xx = np.indices((512, 512))
    for i in range(5):
        # Synthetic horizontal retinal bands, speckle, and dark circular objects.
        rng = np.random.default_rng(20261104 + i)
        image = 170 + 18 * np.cos(yy / 13) + rng.normal(0, 7, (512, 512))
        image[300:390] -= 35
        for x, y, radius in ((65 + i, 335, 13), (170, 354 + i, 18),
                             (300 + i, 331, 20), (405, 361 - i, 14)):
            image[(xx - x) ** 2 + (yy - y) ** 2 <= radius ** 2] -= 95
        Image.fromarray(np.clip(image, 0, 255).astype(np.uint8), "L").save(
            cube / f"{i:04d}.jpg", quality=95
        )
        (seg / f"{i:04d}.svg").write_text(SVG, encoding="utf-8")
    (cube / "info.json").write_text(
        json.dumps({"spacing": [1.0, 1.0, 0.01, 0.01]}, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    generate_example()
    print(f"Generated synthetic example in {EXAMPLE_DIR}")
