# Fourier Coherence Score

**Fourier Coherence Score (FCS)** is an ImageJ/Fiji plugin for quantitative
image-quality assessment. The score ranges from **0 to 1**, with higher values
indicating better image quality.

The plugin supports single images and image stacks. It can generate spatially
resolved image-quality maps, calculate FCS across a stack, identify the slice
with the highest FCS, map regional focus depth, and reconstruct an all-in-focus
image.

## Features

### Single image

- Calculates FCS in overlapping image regions.
- Generates a spatially resolved image-quality map.
- Supports region sizes from 64 x 64 to 1024 x 1024 pixels.
- Provides optional regularization and an optional regional results table.

### Image stack: FCS across stack

- Calculates one FCS value for every slice.
- Displays the FCS-versus-slice plot and results table.
- Reports and displays the slice with the highest FCS.

### Image stack: region-wise refocusing

- Calculates local FCS throughout every slice.
- Generates a spatially resolved focus-depth map.
- Generates a regional FCS image-quality map.
- Reconstructs an all-in-focus image from the stack.
- Provides an optional regional results table.

## Requirements

- ImageJ 1.x or Fiji.
- A grayscale image or grayscale image stack.
- JTransforms 2.4 using the `edu.emory.mathcs.jtransforms` namespace.

In Fiji, the dependency is commonly installed as:

```text
Fiji.app/jars/jtransforms.jar
```

> **Compatibility note:** `JTransforms-3.1-with-dependencies.jar` normally uses
> the newer `org.jtransforms` namespace and is not directly compatible with the
> current plugin source.

## Installation in Fiji

1. Download [Fourier-Coherence-Score.jar](./Fourier-Coherence-Score.jar).
2. Close Fiji.
3. Copy the JAR into the Fiji plugins folder:

   ```text
   Fiji.app/plugins/
   ```

4. Confirm that `jtransforms.jar` is present in `Fiji.app/jars/`.
5. Restart Fiji.
6. Open a grayscale image or image stack.
7. Select **Plugins > Fourier Coherence Score**.

The plugin uses ImageJ/Fiji's built-in **Fire** lookup table. A separate
`jet.lut` file is not required.

## Usage

### Analyze a single image

1. Open a grayscale image.
2. Select **Plugins > Fourier Coherence Score**.
3. Choose the region size and largest feature included.
4. Enable the regularization term if needed.
5. Select whether to show the regional results table.
6. Select **OK**.

The plugin generates an `FCS_Image_Quality_Map` image. Higher map values
indicate regions with better image quality.

### Analyze an image stack

1. Open a grayscale image stack.
2. Select **Plugins > Fourier Coherence Score**.
3. Choose one of the following modes:

   - **FCS across stack**: calculates one score per slice and displays the
     highest-FCS slice.
   - **Region-wise refocusing**: generates the focus-depth map, regional image-
     quality map, and all-in-focus reconstruction.

4. Configure the displayed parameters and select **OK**.

## Main parameters

| Parameter | Description |
| --- | --- |
| Region size | Width and height of each local analysis window. Smaller regions provide finer spatial localization; larger regions provide more stable Fourier statistics. |
| Largest feature included | Largest spatial feature included in the evaluated Fourier-frequency band, expressed in pixels. |
| Enable regularization term | Applies virtual-noise regularization to reduce artificially high scores in weak-signal regions. |
| Show regional results table | Displays the regional measurements used to generate the map. |

Regional windows use 50% overlap. The regional results are interpolated and
smoothed to generate maps with the same width and height as the input image.

## Output names

Depending on the selected mode, the plugin produces:

- `FCS_Image_Quality_Map`
- `FCS_Regional_Results`
- `FCS Across Stack Results`
- `Highest FCS Slice`
- `Focus_Depth_Map`
- `FCS_AllInFocus`

## Source code

The editable plugin source is available in
[Fourier_Coherence_Score_.java](./Fourier_Coherence_Score_.java).

## Author

Hongqiang Ma

## License

This project is released under the [MIT License](./LICENSE).
