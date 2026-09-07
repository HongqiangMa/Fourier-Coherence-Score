/**
 * FCS
 *
 * Standard ImageJ 1.x plugin for full-view stack analysis, regional
 * image-quality mapping, and region-wise refocusing using the
 * Fourier Coherence Score (FCS).
 *
 *
 * Single-image outputs:
 * 1. Smoothed regional FCS score map.
 * 2. Region-wise ResultsTable.
 *
 * Full-view stack outputs:
 * 1. FCS-by-slice ResultsTable and plot.
 * 2. Image copied from the slice with the highest FCS.
 *
 * Region-wise stack outputs:
 * 1. All-in-focus image reconstructed with continuous three-slice axial
 *    weighting.
 * 2. Continuous focus-position map in slice-index units.
 * 3. Smoothed maximum-FCS quality map.
 * 4. Optional region-wise ResultsTable.
 *
 * Reliable regions retain their measured best-focus slices. Unreliable
 * regions receive the focal depth of the nearest reliable region before the
 * continuous focus map and all-in-focus image are reconstructed.
 * A region is reliable only when its maximum FCS is strictly greater than the
 * selected reliability threshold.
 *
 * Requirements:
 * - Standard ImageJ 1.x
 * - JTransforms (jtransforms-*.jar in the ImageJ jars or plugins folder)
 *
 * @author Hongqiang Ma
 * @version 5.8 (September 7, 2026)
 *
 * Copyright (c) 2026 Hongqiang Ma
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextField;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jtransforms.fft.FloatFFT_2D;

public class Fourier_Coherence_Score_ implements PlugIn {

    private static final int DEFAULT_FEATURE_SIZE = 16;
    private static final String[] SINGLE_IMAGE_REGION_SIZE_LABELS = {
            "64 x 64",
            "128 x 128",
            "256 x 256",
            "512 x 512",
            "1024 x 1024"
    };
    private static final int[] SINGLE_IMAGE_REGION_SIZES = {
            64, 128, 256, 512, 1024
    };
    private static final String DEFAULT_SINGLE_IMAGE_REGION_SIZE = "256 x 256";
    private static final double DEFAULT_VIRTUAL_NOISE_STRENGTH = 1.0;
    private static final double DEFAULT_RELIABILITY_THRESHOLD = 0.05;
    private static final double DEFAULT_MAP_SIGMA_FRACTION = 0.5;
    private static final double DEFAULT_AXIAL_SIGMA = 0.5;
    private static final double GAUSSIAN_ACCURACY = 0.01;
    private static final String[] RECONSTRUCTION_METHODS = {
            "Continuous 3-slice weighted",
            "Nearest integer slice"
    };
    private static final String[] STACK_ANALYSIS_MODES = {
            "FCS across stack",
            "Region-wise refocusing"
    };
    private static final int FULL_VIEW_MODE = 0;
    private static final int REGIONAL_MODE = 1;
    private static final double EPSILON = 1.0e-30;

    @Override
    public void run(String arg) {
        ImagePlus input = WindowManager.getCurrentImage();
        if (input == null) {
            IJ.noImage();
            return;
        }



        if (input.getBitDepth() == 24) {
            IJ.error(
                    "Fourier Coherence Score",
                    "RGB images are not supported. Use a grayscale image or z-stack."
            );
            return;
        }

        if (input.getWidth() < 4 || input.getHeight() < 4) {
            IJ.error("Fourier Coherence Score", "The input image must be at least 4 x 4 pixels.");
            return;
        }

        if (input.getStackSize() == 1) {
            runRegionalAnalysis(input);
            return;
        }

        GenericDialog modeDialog = new GenericDialog("Fourier Coherence Score");
        modeDialog.addMessage(
                "Choose how to analyze the current image stack.\n\n"
                        + "Full-view mode calculates one FCS value for each slice.\n"
                        + "Regional mode maps local image quality and focus depth."
        );
        Choice stackModeChoice = new Choice();
        for (String stackAnalysisMode : STACK_ANALYSIS_MODES) {
            stackModeChoice.add(stackAnalysisMode);
        }
        stackModeChoice.select(STACK_ANALYSIS_MODES[FULL_VIEW_MODE]);
        Panel stackModePanel = new Panel(new GridBagLayout());
        addAlignedControlRow(
                stackModePanel,
                0,
                "Analysis mode:",
                stackModeChoice
        );
        addAlignedPanel(modeDialog, stackModePanel);
        modeDialog.showDialog();
        if (modeDialog.wasCanceled()) {
            return;
        }

        int analysisMode = stackModeChoice.getSelectedIndex();
        if (analysisMode == FULL_VIEW_MODE) {
            processFullViewStack(input);
        } else {
            runRegionalAnalysis(input);
        }
    }

    private void runRegionalAnalysis(ImagePlus input) {
        final boolean singleImage = input.getStackSize() == 1;

        GenericDialog dialog = new GenericDialog(
                "Fourier Coherence Score"
        );
        dialog.addMessage(
                singleImage
                        ? "Generate a spatially resolved image-quality map\n"
                        + "for the current image.\n\n"
                        + "The Fourier Coherence Score (FCS) ranges from 0 to 1.\n"
                        + "Higher values indicate better image quality."
                        : "Generate spatially resolved focus-depth maps\n"
                        + "for the current image stack, together with\n"
                        + "an all-in-focus image."
        );

        SingleImageDialogControls singleControls = null;
        RegionalStackDialogControls stackControls = null;
        if (singleImage) {
            singleControls = addSingleImageControls(dialog);
        } else {
            stackControls = addRegionalStackControls(dialog);
        }
        dialog.showDialog();
        if (dialog.wasCanceled()) {
            return;
        }

        int regionWidth;
        int regionHeight;
        int featureSize;
        if (singleImage) {
            int regionSize = SINGLE_IMAGE_REGION_SIZES[
                    singleControls.regionSizeChoice.getSelectedIndex()
            ];
            regionWidth = regionSize;
            regionHeight = regionSize;
            try {
                featureSize = (int) Math.round(Double.parseDouble(
                        singleControls.featureSizeField.getText().trim()
                ));
            } catch (NumberFormatException exception) {
                IJ.error(
                        "Fourier Coherence Score",
                        "Largest feature included must be a valid number."
                );
                return;
            }
        } else {
            int regionSize = SINGLE_IMAGE_REGION_SIZES[
                    stackControls.regionSizeChoice.getSelectedIndex()
            ];
            regionWidth = regionSize;
            regionHeight = regionSize;
            try {
                featureSize = (int) Math.round(Double.parseDouble(
                        stackControls.featureSizeField.getText().trim()
                ));
            } catch (NumberFormatException exception) {
                IJ.error(
                        "Fourier Coherence Score",
                        "Largest feature included must be a valid number."
                );
                return;
            }
        }
        boolean enableVirtualNoise = singleImage
                ? singleControls.regularizationCheckbox.getState()
                : stackControls.regularizationCheckbox.getState();
        double virtualNoiseStrength = DEFAULT_VIRTUAL_NOISE_STRENGTH;
        double reliabilityThreshold = DEFAULT_RELIABILITY_THRESHOLD;
        double mapSigmaFraction = DEFAULT_MAP_SIGMA_FRACTION;
        boolean showResultsTable = singleImage
                ? singleControls.resultsTableCheckbox.getState()
                : stackControls.resultsTableCheckbox.getState();
        boolean writeLog = false;
        int reconstructionMethod = 0;
        double axialSigma = DEFAULT_AXIAL_SIGMA;

        if (regionWidth < 4 || regionHeight < 4) {
            IJ.error(
                    "Fourier Coherence Score",
                    "Region width and height must each be at least 4 pixels."
            );
            return;
        }

        if (featureSize < 4) {
            IJ.error("Fourier Coherence Score", "Feature size must be at least 4 pixels.");
            return;
        }

        if (virtualNoiseStrength < 0.0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "Virtual noise strength must be zero or greater."
            );
            return;
        }

        if (reliabilityThreshold < 0.0 || reliabilityThreshold > 1.0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "The reliability threshold must be between 0 and 1."
            );
            return;
        }

        if (mapSigmaFraction <= 0.0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "The map-smoothing sigma fraction must be greater than zero."
            );
            return;
        }

        if (axialSigma <= 0.0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "The axial-blending sigma must be greater than zero."
            );
            return;
        }

        if (singleImage) {
            processSingleImage(
                    input,
                    regionWidth,
                    regionHeight,
                    featureSize,
                    enableVirtualNoise,
                    virtualNoiseStrength,
                    reliabilityThreshold,
                    mapSigmaFraction,
                    showResultsTable,
                    writeLog
            );
        } else {
            processStack(
                    input,
                    regionWidth,
                    regionHeight,
                    featureSize,
                    enableVirtualNoise,
                    virtualNoiseStrength,
                    reliabilityThreshold,
                    mapSigmaFraction,
                    reconstructionMethod,
                    axialSigma,
                    showResultsTable,
                    writeLog
            );
        }
    }

    private SingleImageDialogControls addSingleImageControls(GenericDialog dialog) {
        Panel controlsPanel = new Panel(new GridBagLayout());
        Choice regionSizeChoice = createRegionSizeChoice();

        TextField featureSizeField = new TextField(
                Integer.toString(DEFAULT_FEATURE_SIZE),
                8
        );
        Checkbox regularizationCheckbox = new Checkbox("", false);
        Checkbox resultsTableCheckbox = new Checkbox("", false);

        addAlignedControlRow(controlsPanel, 0, "Region size (pixels):", regionSizeChoice);
        addAlignedControlRow(
                controlsPanel,
                1,
                "Largest feature included (pixels):",
                featureSizeField
        );
        addAlignedControlRow(
                controlsPanel,
                2,
                "Enable regularization term:",
                regularizationCheckbox
        );
        addAlignedControlRow(
                controlsPanel,
                3,
                "Show regional results table:",
                resultsTableCheckbox
        );

        addAlignedPanel(dialog, controlsPanel);
        return new SingleImageDialogControls(
                regionSizeChoice,
                featureSizeField,
                regularizationCheckbox,
                resultsTableCheckbox
        );
    }

    private RegionalStackDialogControls addRegionalStackControls(GenericDialog dialog) {
        Panel controlsPanel = new Panel(new GridBagLayout());
        Choice regionSizeChoice = createRegionSizeChoice();
        TextField featureSizeField = new TextField(
                Integer.toString(DEFAULT_FEATURE_SIZE),
                8
        );
        Checkbox regularizationCheckbox = new Checkbox("", true);
        Checkbox resultsTableCheckbox = new Checkbox("", false);

        addAlignedControlRow(controlsPanel, 0, "Region size (pixels):", regionSizeChoice);
        addAlignedControlRow(
                controlsPanel,
                1,
                "Largest feature included (pixels):",
                featureSizeField
        );
        addAlignedControlRow(
                controlsPanel,
                2,
                "Enable regularization term:",
                regularizationCheckbox
        );
        addAlignedControlRow(
                controlsPanel,
                3,
                "Show regional results table:",
                resultsTableCheckbox
        );
        addAlignedPanel(dialog, controlsPanel);

        return new RegionalStackDialogControls(
                regionSizeChoice,
                featureSizeField,
                regularizationCheckbox,
                resultsTableCheckbox
        );
    }

    private Choice createRegionSizeChoice() {
        Choice regionSizeChoice = new Choice();
        for (String regionSizeLabel : SINGLE_IMAGE_REGION_SIZE_LABELS) {
            regionSizeChoice.add(regionSizeLabel);
        }
        regionSizeChoice.select(DEFAULT_SINGLE_IMAGE_REGION_SIZE);
        return regionSizeChoice;
    }

    private void addAlignedPanel(GenericDialog dialog, Panel panel) {
        dialog.addPanel(
                panel,
                GridBagConstraints.EAST,
                new Insets(8, 0, 0, 0)
        );
    }

    private void addAlignedControlRow(
            Panel panel,
            int row,
            String labelText,
            java.awt.Component control
    ) {
        GridBagConstraints labelConstraints = new GridBagConstraints();
        labelConstraints.gridx = 0;
        labelConstraints.gridy = row;
        labelConstraints.anchor = GridBagConstraints.EAST;
        labelConstraints.insets = new Insets(4, 4, 4, 16);

        GridBagConstraints controlConstraints = new GridBagConstraints();
        controlConstraints.gridx = 1;
        controlConstraints.gridy = row;
        controlConstraints.anchor = GridBagConstraints.EAST;
        controlConstraints.insets = new Insets(4, 4, 4, 4);

        panel.add(new Label(labelText, Label.RIGHT), labelConstraints);
        panel.add(control, controlConstraints);
    }

    private FullViewDialogControls addFullViewStackControls(GenericDialog dialog) {
        Panel controlsPanel = new Panel(new GridBagLayout());
        TextField featureSizeField = new TextField(
                Integer.toString(DEFAULT_FEATURE_SIZE),
                8
        );
        Checkbox regularizationCheckbox = new Checkbox("", true);

        addAlignedControlRow(
                controlsPanel,
                0,
                "Largest feature included (pixels):",
                featureSizeField
        );
        addAlignedControlRow(
                controlsPanel,
                1,
                "Enable regularization term:",
                regularizationCheckbox
        );
        addAlignedPanel(dialog, controlsPanel);

        return new FullViewDialogControls(
                featureSizeField,
                regularizationCheckbox
        );
    }

    private void processFullViewStack(ImagePlus input) {
        GenericDialog dialog = new GenericDialog("Fourier Coherence Score");
        dialog.addMessage(
                "Calculate the Fourier Coherence Score (FCS)\n"
                        + "for each slice of the current image stack.\n\n"
                        + "The FCS ranges from 0 to 1, with higher values\n"
                        + "indicating better image quality. The slice with\n"
                        + "the highest FCS is reported and displayed."
        );
        FullViewDialogControls controls = addFullViewStackControls(dialog);
        dialog.showDialog();
        if (dialog.wasCanceled()) {
            return;
        }

        int featureSize;
        try {
            featureSize = (int) Math.round(Double.parseDouble(
                    controls.featureSizeField.getText().trim()
            ));
        } catch (NumberFormatException exception) {
            IJ.error(
                    "Fourier Coherence Score",
                    "Largest feature included must be a valid number."
            );
            return;
        }
        double tukeyAlpha = 0.25;
        boolean enableVirtualNoise = controls.regularizationCheckbox.getState();
        boolean writeLog = false;
        double virtualNoiseStrength = DEFAULT_VIRTUAL_NOISE_STRENGTH;

        if (featureSize < 4) {
            IJ.error("Fourier Coherence Score", "Feature size must be at least 4 pixels.");
            return;
        }
        if (!Double.isFinite(tukeyAlpha) || tukeyAlpha < 0.0 || tukeyAlpha > 1.0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "Tukey window alpha must be between 0 and 1."
            );
            return;
        }
        if (!Double.isFinite(virtualNoiseStrength) || virtualNoiseStrength < 0.0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "Virtual noise strength must be zero or greater."
            );
            return;
        }

        int analysisWidth = input.getWidth() - input.getWidth() % 2;
        int analysisHeight = input.getHeight() - input.getHeight() % 2;
        if (analysisWidth < 4 || analysisHeight < 4) {
            IJ.error(
                    "Fourier Coherence Score",
                    "The analysis area must be at least 4 x 4 pixels."
            );
            return;
        }

        int splitWidth = analysisWidth / 2;
        int splitHeight = analysisHeight / 2;
        int sliceCount = input.getStackSize();
        FloatFFT_2D fft = new FloatFFT_2D(splitHeight, splitWidth);
        boolean[] frequencyMask = createFullViewFrequencyMask(
                splitWidth,
                splitHeight,
                featureSize
        );
        int frequencyCount = countTrue(frequencyMask);
        if (frequencyCount == 0) {
            IJ.error(
                    "Fourier Coherence Score",
                    "The selected feature size leaves no Fourier samples in the band."
            );
            return;
        }

        double[] window = createTukeyWindow(splitWidth, splitHeight, tukeyAlpha);
        SpectralMeasurement[] measurements = new SpectralMeasurement[sliceCount];
        ImageStack stack = input.getStack();

        for (int slice = 1; slice <= sliceCount; slice++) {
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                IJ.showStatus("Fourier Coherence Score analysis canceled");
                return;
            }
            measurements[slice - 1] = measureFullViewSpectrum(
                    stack.getProcessor(slice),
                    analysisWidth,
                    analysisHeight,
                    fft,
                    frequencyMask,
                    window
            );
            IJ.showProgress(slice, sliceCount);
            IJ.showStatus(
                    "Measuring full-view Fourier Coherence Score: "
                            + slice + "/" + sliceCount
            );
        }

        double referenceBandPower = 0.0;
        for (SpectralMeasurement measurement : measurements) {
            referenceBandPower = Math.max(
                    referenceBandPower,
                    measurement.geometricBandPower()
            );
        }
        double virtualNoisePower = enableVirtualNoise
                ? virtualNoiseStrength * referenceBandPower
                : 0.0;

        double[] scores = new double[sliceCount];
        double[] sliceAxis = new double[sliceCount];
        int peakSlice = 1;
        double peakScore = -1.0;
        ResultsTable results = new ResultsTable();

        for (int index = 0; index < sliceCount; index++) {
            scores[index] = calculateModifiedFcs(measurements[index], virtualNoisePower);
            sliceAxis[index] = index + 1;
            if (scores[index] > peakScore) {
                peakScore = scores[index];
                peakSlice = index + 1;
            }
            results.incrementCounter();
            results.addValue("Slice", index + 1);
            results.addValue("FCS", scores[index]);
        }

        String prefix = input.getShortTitle();
        results.show(WindowManager.makeUniqueName(prefix + " - FCS Across Stack Results"));
        showFullViewPlot(prefix, sliceAxis, scores, enableVirtualNoise, peakSlice, peakScore);

        ImagePlus highestFcsSlice = new ImagePlus(
                WindowManager.makeUniqueName(
                        prefix + " - Highest FCS Slice " + peakSlice
                ),
                stack.getProcessor(peakSlice).duplicate()
        );
        highestFcsSlice.setCalibration(input.getCalibration().copy());
        highestFcsSlice.setDisplayRange(
                input.getDisplayRangeMin(),
                input.getDisplayRangeMax()
        );
        highestFcsSlice.show();

        if (writeLog) {
            IJ.log("Fourier Coherence Score - full-view stack analysis");
            IJ.log("  Image: " + input.getTitle());
            IJ.log("  Analysis area: 0, 0, " + analysisWidth + " x " + analysisHeight + " px");
            IJ.log("  Largest feature included: " + featureSize + " px");
            IJ.log("  Fourier samples in band: " + frequencyCount);
            IJ.log("  Tukey alpha: " + IJ.d2s(tukeyAlpha, 3));
            IJ.log("  Virtual noise regularization: "
                    + (enableVirtualNoise ? "enabled" : "disabled"));
            if (enableVirtualNoise) {
                IJ.log("  Virtual noise strength: " + IJ.d2s(virtualNoiseStrength, 4));
                IJ.log("  Reference band power: " + IJ.d2s(referenceBandPower, 6));
                IJ.log("  Virtual noise power: " + IJ.d2s(virtualNoisePower, 6));
            }
            IJ.log(
                    "  Highest-FCS slice: " + peakSlice
                            + ", FCS=" + IJ.d2s(peakScore, 6)
            );
        }

        IJ.showProgress(1.0);
        IJ.showStatus(
                "FCS across stack complete; highest FCS at slice " + peakSlice
        );
    }

    private boolean[] createFullViewFrequencyMask(
            int splitWidth,
            int splitHeight,
            int featureSize
    ) {
        boolean[] mask = new boolean[splitWidth * splitHeight];
        double cutoffCyclesPerInputPixel = 1.0 / featureSize;

        for (int y = 0; y < splitHeight; y++) {
            int ky = y <= splitHeight / 2 ? y : y - splitHeight;
            double fy = ky / (2.0 * splitHeight);
            for (int x = 0; x < splitWidth; x++) {
                int kx = x <= splitWidth / 2 ? x : x - splitWidth;
                double fx = kx / (2.0 * splitWidth);
                double radialFrequency = Math.sqrt(fx * fx + fy * fy);
                mask[y * splitWidth + x] = radialFrequency >= cutoffCyclesPerInputPixel;
            }
        }
        return mask;
    }

    private SpectralMeasurement measureFullViewSpectrum(
            ImageProcessor processor,
            int analysisWidth,
            int analysisHeight,
            FloatFFT_2D fft,
            boolean[] frequencyMask,
            double[] window
    ) {
        int splitWidth = analysisWidth / 2;
        int splitHeight = analysisHeight / 2;
        int sampleCount = splitWidth * splitHeight;
        float[] imageA = new float[2 * sampleCount];
        float[] imageB = new float[2 * sampleCount];
        double meanA = 0.0;
        double meanB = 0.0;

        for (int y = 0; y < splitHeight; y++) {
            int sourceY = 2 * y;
            for (int x = 0; x < splitWidth; x++) {
                int sourceX = 2 * x;
                float a = 0.5f * (
                        processor.getf(sourceX, sourceY)
                                + processor.getf(sourceX + 1, sourceY + 1)
                );
                float b = 0.5f * (
                        processor.getf(sourceX + 1, sourceY)
                                + processor.getf(sourceX, sourceY + 1)
                );
                int index = 2 * (y * splitWidth + x);
                imageA[index] = a;
                imageB[index] = b;
                meanA += a;
                meanB += b;
            }
        }

        meanA /= sampleCount;
        meanB /= sampleCount;
        for (int index = 0; index < sampleCount; index++) {
            int complexIndex = 2 * index;
            imageA[complexIndex] = (float) ((imageA[complexIndex] - meanA) * window[index]);
            imageB[complexIndex] = (float) ((imageB[complexIndex] - meanB) * window[index]);
        }

        fft.complexForward(imageA);
        fft.complexForward(imageB);

        double crossReal = 0.0;
        double crossImaginary = 0.0;
        double powerA = 0.0;
        double powerB = 0.0;
        for (int index = 0; index < sampleCount; index++) {
            if (!frequencyMask[index]) {
                continue;
            }
            int complexIndex = 2 * index;
            double aReal = imageA[complexIndex];
            double aImaginary = imageA[complexIndex + 1];
            double bReal = imageB[complexIndex];
            double bImaginary = imageB[complexIndex + 1];
            crossReal += aReal * bReal + aImaginary * bImaginary;
            crossImaginary += aImaginary * bReal - aReal * bImaginary;
            powerA += aReal * aReal + aImaginary * aImaginary;
            powerB += bReal * bReal + bImaginary * bImaginary;
        }

        return new SpectralMeasurement(
                Math.hypot(crossReal, crossImaginary),
                powerA,
                powerB
        );
    }

    private double[] createTukeyWindow(int width, int height, double alpha) {
        double[] windowX = createTukeyVector(width, alpha);
        double[] windowY = createTukeyVector(height, alpha);
        double[] window = new double[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                window[y * width + x] = windowX[x] * windowY[y];
            }
        }
        return window;
    }

    private double[] createTukeyVector(int length, double alpha) {
        double[] values = new double[length];
        if (length == 1 || alpha <= 0.0) {
            for (int index = 0; index < length; index++) {
                values[index] = 1.0;
            }
            return values;
        }

        for (int index = 0; index < length; index++) {
            double position = index / (double) (length - 1);
            if (position < alpha / 2.0) {
                values[index] = 0.5 * (1.0 + Math.cos(
                        Math.PI * (2.0 * position / alpha - 1.0)
                ));
            } else if (position <= 1.0 - alpha / 2.0) {
                values[index] = 1.0;
            } else {
                values[index] = 0.5 * (1.0 + Math.cos(
                        Math.PI * (2.0 * position / alpha - 2.0 / alpha + 1.0)
                ));
            }
        }
        return values;
    }

    private int countTrue(boolean[] values) {
        int count = 0;
        for (boolean value : values) {
            if (value) {
                count++;
            }
        }
        return count;
    }

    private void showFullViewPlot(
            String prefix,
            double[] sliceAxis,
            double[] scores,
            boolean regularized,
            int peakSlice,
        double peakScore
    ) {
        String title = WindowManager.makeUniqueName(
                prefix + " - FCS Across Stack"
                        + (regularized ? " (Regularized)" : "")
        );
        Plot plot = new Plot(title, "Slice", "Fourier Coherence Score (FCS)");
        double maximum = 0.0;
        for (double score : scores) {
            maximum = Math.max(maximum, score);
        }
        double upperLimit = maximum > 0.0 ? Math.min(1.05, maximum * 1.1) : 1.0;
        double lastSlice = scores.length > 1 ? scores.length : 2.0;
        plot.setLimits(1.0, lastSlice, 0.0, upperLimit);
        plot.setLineWidth(2);
        plot.setColor("black");
        plot.addPoints(sliceAxis, scores, Plot.LINE);
        plot.setLineWidth(3);
        plot.addPoints(sliceAxis, scores, Plot.CIRCLE);
        plot.setColor("red");
        plot.setLineWidth(5);
        plot.addPoints(
                new double[]{peakSlice},
                new double[]{peakScore},
                Plot.CIRCLE
        );
        plot.show();
    }

    private void processSingleImage(
            ImagePlus input,
            int requestedRegionWidth,
            int requestedRegionHeight,
            int featureSize,
            boolean enableVirtualNoise,
            double virtualNoiseStrength,
            double reliabilityThreshold,
            double mapSigmaFraction,
            boolean showResultsTable,
            boolean writeLog
    ) {
        final int width = input.getWidth();
        final int height = input.getHeight();
        final int regionWidth = Math.min(requestedRegionWidth, width);
        final int regionHeight = Math.min(requestedRegionHeight, height);

        int[] xPositions = createPositions(width, regionWidth);
        int[] yPositions = createPositions(height, regionHeight);
        List<RegionData> regions = createRegions(
                xPositions,
                yPositions,
                regionWidth,
                regionHeight
        );

        Map<Long, FftPlan> planCache = new HashMap<Long, FftPlan>();
        measureSingleImageRegions(
                input.getProcessor(),
                regions,
                featureSize,
                enableVirtualNoise,
                virtualNoiseStrength,
                reliabilityThreshold,
                planCache
        );

        double sigmaX = Math.max(0.5, regionWidth * mapSigmaFraction);
        double sigmaY = Math.max(0.5, regionHeight * mapSigmaFraction);
        FloatProcessor qualityMap = createQualityMap(
                regions,
                width,
                height,
                sigmaX,
                sigmaY
        );
        ResultsTable results = showResultsTable
                ? createSingleImageResultsTable(regions)
                : null;
        showSingleImageOutputs(input, qualityMap, results);

        if (writeLog) {
            IJ.log(
                    "Fourier Coherence Score single-image regional mode: "
                            + regions.size() + " regions, sigma = "
                            + IJ.d2s(sigmaX, 1) + " x " + IJ.d2s(sigmaY, 1) + " px, "
                            + "virtual noise = " + (enableVirtualNoise
                            ? "on (strength " + IJ.d2s(virtualNoiseStrength, 3) + ")"
                            : "off") + ". The map reports relative local FCS quality, "
                            + "not absolute focal depth."
            );
        }
        IJ.showProgress(1.0);
        IJ.showStatus("FCS single-image score map complete");
    }

    private void measureSingleImageRegions(
            ImageProcessor processor,
            List<RegionData> regions,
            int featureSize,
            boolean enableVirtualNoise,
            double virtualNoiseStrength,
            double reliabilityThreshold,
            Map<Long, FftPlan> planCache
    ) {
        SpectralMeasurement[] measurements = new SpectralMeasurement[regions.size()];
        double imageReferenceBandPower = 0.0;

        for (int regionIndex = 0; regionIndex < regions.size(); regionIndex++) {
            RegionData region = regions.get(regionIndex);
            int analysisWidth = region.width - (region.width % 2);
            int analysisHeight = region.height - (region.height % 2);
            FftPlan plan = getPlan(planCache, analysisWidth, analysisHeight, featureSize);
            SpectralMeasurement measurement = measureSpectrum(
                    processor,
                    region.x0,
                    region.y0,
                    plan
            );
            measurements[regionIndex] = measurement;
            imageReferenceBandPower = Math.max(
                    imageReferenceBandPower,
                    measurement.geometricBandPower()
            );
        }

        double virtualNoisePower = enableVirtualNoise
                ? virtualNoiseStrength * imageReferenceBandPower
                : 0.0;

        for (int regionIndex = 0; regionIndex < regions.size(); regionIndex++) {
            RegionData region = regions.get(regionIndex);
            double score = calculateModifiedFcs(
                    measurements[regionIndex],
                    virtualNoisePower
            );
            region.bestSlice = 1;
            region.maximumFcs = score;
            region.reliable = score > reliabilityThreshold;
            region.focusSlice = Double.NaN;
            region.reconstructionSlice = 0;

            IJ.showStatus(
                    "Measuring single-image regional FCS: region " + (regionIndex + 1)
                            + " of " + regions.size()
            );
            IJ.showProgress(regionIndex + 1, regions.size());
        }
    }

    private FloatProcessor createQualityMap(
            List<RegionData> regions,
            int width,
            int height,
            double sigmaX,
            double sigmaY
    ) {
        FloatProcessor qualityNumerator = new FloatProcessor(width, height);
        FloatProcessor allRegionWeight = new FloatProcessor(width, height);

        for (RegionData region : regions) {
            int x = clamp((int) Math.round(region.centerX), 0, width - 1);
            int y = clamp((int) Math.round(region.centerY), 0, height - 1);
            qualityNumerator.setf(
                    x,
                    y,
                    qualityNumerator.getf(x, y) + (float) region.maximumFcs
            );
            allRegionWeight.setf(x, y, allRegionWeight.getf(x, y) + 1.0f);
        }

        GaussianBlur gaussianBlur = new GaussianBlur();
        gaussianBlur.blurGaussian(
                qualityNumerator,
                sigmaX,
                sigmaY,
                GAUSSIAN_ACCURACY
        );
        gaussianBlur.blurGaussian(
                allRegionWeight,
                sigmaX,
                sigmaY,
                GAUSSIAN_ACCURACY
        );

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float weight = allRegionWeight.getf(x, y);
                float score = weight > EPSILON
                        ? qualityNumerator.getf(x, y) / weight
                        : 0.0f;
                qualityNumerator.setf(
                        x,
                        y,
                        (float) Math.max(0.0, Math.min(1.0, score))
                );
            }
        }
        return qualityNumerator;
    }

    private ResultsTable createSingleImageResultsTable(List<RegionData> regions) {
        ResultsTable results = new ResultsTable();
        for (RegionData region : regions) {
            results.incrementCounter();
            results.addValue("Region", region.number);
            results.addValue("Grid_Row", region.gridRow);
            results.addValue("Grid_Column", region.gridColumn);
            results.addValue("X_px", region.x0);
            results.addValue("Y_px", region.y0);
            results.addValue("Width_px", region.width);
            results.addValue("Height_px", region.height);
            results.addValue("FCS", region.maximumFcs);
            results.addValue("Above_Threshold", region.reliable ? 1 : 0);
        }
        return results;
    }

    private void showSingleImageOutputs(
            ImagePlus input,
            FloatProcessor qualityMap,
            ResultsTable results
    ) {
        String prefix = input.getShortTitle();
        Calibration spatialCalibration = input.getCalibration().copy();
        ImagePlus qualityMapImage = new ImagePlus(
                WindowManager.makeUniqueName(prefix + "_FCS_Image_Quality_Map"),
                qualityMap
        );
        Calibration qualityCalibration = spatialCalibration.copy();
        qualityCalibration.setValueUnit("FCS");
        qualityMapImage.setCalibration(qualityCalibration);
        qualityMapImage.resetDisplayRange();
        qualityMapImage.show();
        IJ.run(qualityMapImage, "Jet", "");
        if (results != null) {
            results.show(WindowManager.makeUniqueName(prefix + "_FCS_Regional_Results"));
        }
    }

    private void processStack(
            ImagePlus input,
            int requestedRegionWidth,
            int requestedRegionHeight,
            int featureSize,
            boolean enableVirtualNoise,
            double virtualNoiseStrength,
            double reliabilityThreshold,
            double mapSigmaFraction,
            int reconstructionMethod,
            double axialSigma,
            boolean showResultsTable,
            boolean writeLog
    ) {
        final int width = input.getWidth();
        final int height = input.getHeight();
        final int slices = input.getStackSize();
        final int regionWidth = Math.min(requestedRegionWidth, width);
        final int regionHeight = Math.min(requestedRegionHeight, height);

        int[] xPositions = createPositions(width, regionWidth);
        int[] yPositions = createPositions(height, regionHeight);
        List<RegionData> regions = createRegions(
                xPositions,
                yPositions,
                regionWidth,
                regionHeight
        );

        ImageStack stack = input.getStack();
        Map<Long, FftPlan> planCache = new HashMap<Long, FftPlan>();
        double[] globalScoreSum = new double[slices];

        measureRegions(
                stack,
                regions,
                featureSize,
                enableVirtualNoise,
                virtualNoiseStrength,
                reliabilityThreshold,
                planCache,
                globalScoreSum
        );

        int reliableRegions = countReliableRegions(regions);
        int fallbackSlice = findGlobalBestSlice(globalScoreSum);
        assignRegionalFocusDepths(regions, reliableRegions, fallbackSlice);
        double sigmaX = Math.max(0.5, regionWidth * mapSigmaFraction);
        double sigmaY = Math.max(0.5, regionHeight * mapSigmaFraction);

        Reconstruction reconstruction = reconstructOutputs(
                input,
                stack,
                regions,
                width,
                height,
                slices,
                reliableRegions,
                fallbackSlice,
                sigmaX,
                sigmaY,
                reconstructionMethod,
                axialSigma,
                writeLog
        );
        updateRegionFocusValues(regions, reconstruction.focusMap, slices);
        ResultsTable results = showResultsTable ? createResultsTable(regions) : null;
        showOutputs(input, reconstruction, results, slices);

        int unreliableRegions = regions.size() - reliableRegions;
        if (writeLog) {
            IJ.log(
                    "Fourier Coherence Score regional stack mode: "
                            + reliableRegions + " reliable regions, "
                            + unreliableRegions + " nearest-neighbor assigned regions (threshold = "
                            + IJ.d2s(reliabilityThreshold, 3) + ", sigma = "
                            + IJ.d2s(sigmaX, 1) + " x " + IJ.d2s(sigmaY, 1) + " px, "
                            + "virtual noise = " + (enableVirtualNoise
                            ? "on (strength " + IJ.d2s(virtualNoiseStrength, 3) + ")"
                            : "off") + ", "
                            + "reconstruction = " + RECONSTRUCTION_METHODS[reconstructionMethod]
                            + (reconstructionMethod == 0
                            ? ", axial sigma = " + IJ.d2s(axialSigma, 2) + " slices)."
                            : ").")
            );
        }
        IJ.showProgress(1.0);
        IJ.showStatus("Regional Fourier Coherence Score processing complete");
    }

    private int[] createPositions(int length, int regionSize) {
        if (regionSize >= length) {
            return new int[]{0};
        }

        int step = Math.max(1, regionSize / 2);
        List<Integer> positions = new ArrayList<Integer>();
        int position = 0;

        while (true) {
            positions.add(position);
            if (position + regionSize >= length) {
                break;
            }

            int nextPosition = position + step;
            if (nextPosition + regionSize >= length) {
                nextPosition = length - regionSize;
            }
            if (nextPosition <= position) {
                break;
            }
            position = nextPosition;
        }

        int[] result = new int[positions.size()];
        for (int i = 0; i < positions.size(); i++) {
            result[i] = positions.get(i);
        }
        return result;
    }

    private List<RegionData> createRegions(
            int[] xPositions,
            int[] yPositions,
            int regionWidth,
            int regionHeight
    ) {
        List<RegionData> regions = new ArrayList<RegionData>();
        int regionNumber = 0;

        for (int row = 0; row < yPositions.length; row++) {
            int y0 = yPositions[row];
            for (int column = 0; column < xPositions.length; column++) {
                int x0 = xPositions[column];
                regionNumber++;

                RegionData region = new RegionData(
                        regionNumber,
                        row + 1,
                        column + 1,
                        x0,
                        y0,
                        regionWidth,
                        regionHeight
                );
                regions.add(region);
            }
        }
        return regions;
    }

    private void measureRegions(
            ImageStack stack,
            List<RegionData> regions,
            int featureSize,
            boolean enableVirtualNoise,
            double virtualNoiseStrength,
            double reliabilityThreshold,
            Map<Long, FftPlan> planCache,
            double[] globalScoreSum
    ) {
        int slices = stack.getSize();
        int totalRegions = regions.size();

        for (int regionIndex = 0; regionIndex < totalRegions; regionIndex++) {
            RegionData region = regions.get(regionIndex);
            int analysisWidth = region.width - (region.width % 2);
            int analysisHeight = region.height - (region.height % 2);
            FftPlan plan = getPlan(planCache, analysisWidth, analysisHeight, featureSize);

            SpectralMeasurement[] measurements = new SpectralMeasurement[slices];
            double referenceBandPower = 0.0;
            for (int slice = 1; slice <= slices; slice++) {
                ImageProcessor processor = stack.getProcessor(slice);
                SpectralMeasurement measurement = measureSpectrum(
                        processor,
                        region.x0,
                        region.y0,
                        plan
                );
                measurements[slice - 1] = measurement;
                referenceBandPower = Math.max(
                        referenceBandPower,
                        measurement.geometricBandPower()
                );
            }

            double virtualNoisePower = enableVirtualNoise
                    ? virtualNoiseStrength * referenceBandPower
                    : 0.0;
            int bestSlice = 1;
            double bestScore = -1.0;

            for (int slice = 1; slice <= slices; slice++) {
                double score = calculateModifiedFcs(
                        measurements[slice - 1],
                        virtualNoisePower
                );
                globalScoreSum[slice - 1] += score;

                if (score > bestScore) {
                    bestScore = score;
                    bestSlice = slice;
                }
            }

            if (!Double.isFinite(bestScore) || bestScore < 0.0) {
                bestScore = 0.0;
            }

            region.bestSlice = bestSlice;
            region.maximumFcs = bestScore;
            region.reliable = bestScore > reliabilityThreshold;

            IJ.showStatus(
                    "Measuring regional FCS: region " + (regionIndex + 1)
                            + " of " + totalRegions
            );
            IJ.showProgress(regionIndex + 1, totalRegions);
        }
    }

    private int countReliableRegions(List<RegionData> regions) {
        int reliableCount = 0;
        for (RegionData region : regions) {
            if (region.reliable) {
                reliableCount++;
            }
        }
        return reliableCount;
    }

    private void assignRegionalFocusDepths(
            List<RegionData> regions,
            int reliableRegions,
            int fallbackSlice
    ) {
        if (reliableRegions == 0) {
            for (RegionData region : regions) {
                region.assignedFocusSlice = fallbackSlice;
                region.focusSourceRegion = 0;
                region.focusSourceDistance = Double.NaN;
            }
            return;
        }

        for (RegionData region : regions) {
            if (region.reliable) {
                region.assignedFocusSlice = region.bestSlice;
                region.focusSourceRegion = region.number;
                region.focusSourceDistance = 0.0;
                continue;
            }

            RegionData nearestReliable = null;
            double nearestDistanceSquared = Double.POSITIVE_INFINITY;
            for (RegionData candidate : regions) {
                if (!candidate.reliable) {
                    continue;
                }
                double dx = candidate.centerX - region.centerX;
                double dy = candidate.centerY - region.centerY;
                double distanceSquared = dx * dx + dy * dy;
                if (distanceSquared < nearestDistanceSquared) {
                    nearestDistanceSquared = distanceSquared;
                    nearestReliable = candidate;
                }
            }

            if (nearestReliable != null) {
                region.assignedFocusSlice = nearestReliable.bestSlice;
                region.focusSourceRegion = nearestReliable.number;
                region.focusSourceDistance = Math.sqrt(nearestDistanceSquared);
            } else {
                region.assignedFocusSlice = fallbackSlice;
                region.focusSourceRegion = 0;
                region.focusSourceDistance = Double.NaN;
            }
        }
    }

    private int findGlobalBestSlice(double[] globalScoreSum) {
        int bestSlice = 1;
        double bestScore = -1.0;
        for (int i = 0; i < globalScoreSum.length; i++) {
            if (globalScoreSum[i] > bestScore) {
                bestScore = globalScoreSum[i];
                bestSlice = i + 1;
            }
        }
        return bestSlice;
    }

    private Reconstruction reconstructOutputs(
            ImagePlus input,
            ImageStack stack,
            List<RegionData> regions,
            int width,
            int height,
            int slices,
            int reliableRegions,
            int fallbackSlice,
            double sigmaX,
            double sigmaY,
            int reconstructionMethod,
            double axialSigma,
            boolean writeLog
    ) {
        FloatProcessor focusNumerator = new FloatProcessor(width, height);
        FloatProcessor focusWeight = new FloatProcessor(width, height);
        FloatProcessor validWeight = new FloatProcessor(width, height);
        FloatProcessor qualityNumerator = new FloatProcessor(width, height);
        FloatProcessor allRegionWeight = new FloatProcessor(width, height);

        seedRegionalMeasurements(
                regions,
                focusNumerator,
                focusWeight,
                validWeight,
                qualityNumerator,
                allRegionWeight
        );

        GaussianBlur focusBlur = new GaussianBlur();
        focusBlur.blurGaussian(
                focusNumerator,
                sigmaX,
                sigmaY,
                GAUSSIAN_ACCURACY
        );
        focusBlur.blurGaussian(focusWeight, sigmaX, sigmaY, GAUSSIAN_ACCURACY);
        focusBlur.blurGaussian(validWeight, sigmaX, sigmaY, GAUSSIAN_ACCURACY);

        if (reliableRegions == 0 && writeLog) {
            IJ.log(
                    "Fourier Coherence Score warning: no region passed the reliability threshold; "
                            + "the globally highest-scoring slice was assigned to all regions."
            );
        }

        GaussianBlur gaussianBlur = new GaussianBlur();
        gaussianBlur.blurGaussian(
                qualityNumerator,
                sigmaX,
                sigmaY,
                GAUSSIAN_ACCURACY
        );
        gaussianBlur.blurGaussian(allRegionWeight, sigmaX, sigmaY, GAUSSIAN_ACCURACY);

        byte[] reliabilityPixels = new byte[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float allWeight = allRegionWeight.getf(x, y);

                float localFocusWeight = focusWeight.getf(x, y);
                float focusValue = localFocusWeight > EPSILON
                        ? focusNumerator.getf(x, y) / localFocusWeight
                        : fallbackSlice;
                focusNumerator.setf(
                        x,
                        y,
                        (float) Math.max(1.0, Math.min(slices, focusValue))
                );

                float qualityValue = allWeight > EPSILON
                        ? qualityNumerator.getf(x, y) / allWeight
                        : 0.0f;
                qualityNumerator.setf(
                        x,
                        y,
                        (float) Math.max(0.0, Math.min(1.0, qualityValue))
                );

                double reliability = allWeight > EPSILON
                        ? validWeight.getf(x, y) / allWeight
                        : 0.0;
                reliability = Math.max(0.0, Math.min(1.0, reliability));
                reliabilityPixels[y * width + x] = (byte) Math.round(255.0 * reliability);
            }
        }

        ByteProcessor reliabilityMap = new ByteProcessor(width, height, reliabilityPixels, null);
        ImageProcessor allInFocus = reconstructFromFocusMap(
                input,
                stack,
                focusNumerator,
                slices,
                reconstructionMethod,
                axialSigma
        );

        return new Reconstruction(
                allInFocus,
                focusNumerator,
                qualityNumerator,
                reliabilityMap
        );
    }

    private void seedRegionalMeasurements(
            List<RegionData> regions,
            FloatProcessor focusNumerator,
            FloatProcessor focusWeight,
            FloatProcessor validWeight,
            FloatProcessor qualityNumerator,
            FloatProcessor allRegionWeight
    ) {
        for (RegionData region : regions) {
            int x = clamp((int) Math.round(region.centerX), 0, focusNumerator.getWidth() - 1);
            int y = clamp((int) Math.round(region.centerY), 0, focusNumerator.getHeight() - 1);

            qualityNumerator.setf(x, y, qualityNumerator.getf(x, y) + (float) region.maximumFcs);
            allRegionWeight.setf(x, y, allRegionWeight.getf(x, y) + 1.0f);

            focusNumerator.setf(
                    x,
                    y,
                    focusNumerator.getf(x, y) + (float) region.assignedFocusSlice
            );
            focusWeight.setf(x, y, focusWeight.getf(x, y) + 1.0f);

            if (region.reliable) {
                validWeight.setf(x, y, validWeight.getf(x, y) + 1.0f);
            }
        }
    }

    private ImageProcessor reconstructFromFocusMap(
            ImagePlus input,
            ImageStack stack,
            FloatProcessor focusMap,
            int slices,
            int reconstructionMethod,
            double axialSigma
    ) {
        int width = input.getWidth();
        int height = input.getHeight();
        ImageProcessor allInFocus = input.getProcessor().createProcessor(width, height);
        ImageProcessor[] processors = new ImageProcessor[slices];
        for (int slice = 1; slice <= slices; slice++) {
            processors[slice - 1] = stack.getProcessor(slice);
        }

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double focusDepth = Math.max(1.0, Math.min(slices, focusMap.getf(x, y)));

                if (reconstructionMethod == 0) {
                    allInFocus.setf(
                            x,
                            y,
                            calculateAxiallyWeightedIntensity(
                                    processors,
                                    x,
                                    y,
                                    focusDepth,
                                    slices,
                                    axialSigma
                            )
                    );
                } else {
                    int slice = clamp((int) Math.round(focusDepth), 1, slices);
                    allInFocus.setf(x, y, processors[slice - 1].getf(x, y));
                }
            }
            IJ.showProgress(y + 1, height);
        }
        return allInFocus;
    }

    private float calculateAxiallyWeightedIntensity(
            ImageProcessor[] processors,
            int x,
            int y,
            double focusDepth,
            int slices,
            double axialSigma
    ) {
        int centerSlice = clamp((int) Math.round(focusDepth), 1, slices);
        int firstSlice = Math.max(1, centerSlice - 1);
        int lastSlice = Math.min(slices, centerSlice + 1);
        double inverseTwoSigmaSquared = 1.0 / (2.0 * axialSigma * axialSigma);
        double weightedIntensity = 0.0;
        double weightSum = 0.0;

        for (int slice = firstSlice; slice <= lastSlice; slice++) {
            double distance = slice - focusDepth;
            double weight = Math.exp(-distance * distance * inverseTwoSigmaSquared);
            weightedIntensity += weight * processors[slice - 1].getf(x, y);
            weightSum += weight;
        }

        if (weightSum <= EPSILON) {
            return processors[centerSlice - 1].getf(x, y);
        }
        return (float) (weightedIntensity / weightSum);
    }

    private void updateRegionFocusValues(
            List<RegionData> regions,
            FloatProcessor focusMap,
            int slices
    ) {
        for (RegionData region : regions) {
            int x = clamp((int) Math.round(region.centerX), 0, focusMap.getWidth() - 1);
            int y = clamp((int) Math.round(region.centerY), 0, focusMap.getHeight() - 1);
            region.focusSlice = Math.max(1.0, Math.min(slices, focusMap.getf(x, y)));
            region.reconstructionSlice = clamp(
                    (int) Math.round(region.focusSlice),
                    1,
                    slices
            );
        }
    }

    private int clamp(int value, int minimum, int maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    private ResultsTable createResultsTable(List<RegionData> regions) {
        ResultsTable results = new ResultsTable();
        for (RegionData region : regions) {
            results.incrementCounter();
            results.addValue("Region", region.number);
            results.addValue("Grid_Row", region.gridRow);
            results.addValue("Grid_Column", region.gridColumn);
            results.addValue("X_px", region.x0);
            results.addValue("Y_px", region.y0);
            results.addValue("Width_px", region.width);
            results.addValue("Height_px", region.height);
            results.addValue("Best_Slice", region.bestSlice);
            results.addValue("Maximum_FCS", region.maximumFcs);
            results.addValue("Reliable", region.reliable ? 1 : 0);
            results.addValue("Assigned_Focus_Slice", region.assignedFocusSlice);
            results.addValue("Focus_Source_Region", region.focusSourceRegion);
            results.addValue("Focus_Source_Distance_px", region.focusSourceDistance);
            results.addValue("Smoothed_Focus_Slice", region.focusSlice);
            results.addValue("Reconstruction_Slice", region.reconstructionSlice);
        }
        return results;
    }

    private FftPlan getPlan(
            Map<Long, FftPlan> cache,
            int analysisWidth,
            int analysisHeight,
            int featureSize
    ) {
        long key = (((long) analysisWidth) << 32) | (analysisHeight & 0xffffffffL);
        FftPlan plan = cache.get(key);
        if (plan == null) {
            plan = new FftPlan(analysisWidth, analysisHeight, featureSize);
            cache.put(key, plan);
        }
        return plan;
    }

    private SpectralMeasurement measureSpectrum(
            ImageProcessor processor,
            int x0,
            int y0,
            FftPlan plan
    ) {
        final int splitWidth = plan.splitWidth;
        final int splitHeight = plan.splitHeight;
        final float[] imageA = plan.imageA;
        final float[] imageB = plan.imageB;

        for (int y = 0; y < splitHeight; y++) {
            int sourceY = y0 + 2 * y;
            for (int x = 0; x < splitWidth; x++) {
                int sourceX = x0 + 2 * x;
                int index = 2 * (y * splitWidth + x);

                float topLeft = processor.getf(sourceX, sourceY);
                float topRight = processor.getf(sourceX + 1, sourceY);
                float bottomLeft = processor.getf(sourceX, sourceY + 1);
                float bottomRight = processor.getf(sourceX + 1, sourceY + 1);

                imageA[index] = 0.5f * (topLeft + bottomRight);
                imageA[index + 1] = 0.0f;
                imageB[index] = 0.5f * (topRight + bottomLeft);
                imageB[index + 1] = 0.0f;
            }
        }

        plan.fft.complexForward(imageA);
        plan.fft.complexForward(imageB);

        double numeratorReal = 0.0;
        double numeratorImaginary = 0.0;
        double denominatorA = 0.0;
        double denominatorB = 0.0;

        for (int bin : plan.includedBins) {
            int index = 2 * bin;
            double aReal = imageA[index];
            double aImaginary = imageA[index + 1];
            double bReal = imageB[index];
            double bImaginary = imageB[index + 1];

            numeratorReal += aReal * bReal + aImaginary * bImaginary;
            numeratorImaginary += aImaginary * bReal - aReal * bImaginary;
            denominatorA += aReal * aReal + aImaginary * aImaginary;
            denominatorB += bReal * bReal + bImaginary * bImaginary;
        }

        double crossMagnitude = Math.hypot(numeratorReal, numeratorImaginary);
        return new SpectralMeasurement(crossMagnitude, denominatorA, denominatorB);
    }

    private double calculateModifiedFcs(
            SpectralMeasurement measurement,
            double virtualNoisePower
    ) {
        double denominator = Math.sqrt(
                (measurement.powerA + virtualNoisePower)
                        * (measurement.powerB + virtualNoisePower)
        );
        if (!Double.isFinite(denominator) || denominator <= EPSILON) {
            return 0.0;
        }

        double score = measurement.crossMagnitude / denominator;
        if (!Double.isFinite(score)) {
            return 0.0;
        }
        return Math.max(0.0, Math.min(1.0, score));
    }

    private void showOutputs(
            ImagePlus input,
            Reconstruction reconstruction,
            ResultsTable results,
            int slices
    ) {
        String prefix = input.getShortTitle();
        Calibration spatialCalibration = input.getCalibration().copy();

        ImagePlus allInFocusImage = new ImagePlus(
                WindowManager.makeUniqueName(prefix + "_FCS_AllInFocus"),
                reconstruction.allInFocus
        );
        allInFocusImage.setCalibration(spatialCalibration.copy());
        allInFocusImage.setDisplayRange(input.getDisplayRangeMin(), input.getDisplayRangeMax());
        allInFocusImage.show();

        ImagePlus focusMapImage = new ImagePlus(
                WindowManager.makeUniqueName(prefix + "_Focus_Depth_Map"),
                reconstruction.focusMap
        );
        Calibration focusCalibration = spatialCalibration.copy();
        focusCalibration.setValueUnit("slice index");
        focusMapImage.setCalibration(focusCalibration);
        focusMapImage.setDisplayRange(1, slices);
        focusMapImage.resetDisplayRange();
        focusMapImage.show();
        IJ.run(focusMapImage, "Jet", "");

        ImagePlus qualityMapImage = new ImagePlus(
                WindowManager.makeUniqueName(prefix + "_FCS_Image_Quality_Map"),
                reconstruction.qualityMap
        );
        Calibration qualityCalibration = spatialCalibration.copy();
        qualityCalibration.setValueUnit("FCS");
        qualityMapImage.setCalibration(qualityCalibration);
        qualityMapImage.setDisplayRange(0.0, 1.0);
        qualityMapImage.resetDisplayRange();
        qualityMapImage.show();
        IJ.run(qualityMapImage, "Jet", "");

 

        if (results != null) {
            results.show(WindowManager.makeUniqueName(prefix + "_FCS_Regional_Results"));
        }
    }

    private static class FullViewDialogControls {
        final TextField featureSizeField;
        final Checkbox regularizationCheckbox;

        FullViewDialogControls(
                TextField featureSizeField,
                Checkbox regularizationCheckbox
        ) {
            this.featureSizeField = featureSizeField;
            this.regularizationCheckbox = regularizationCheckbox;
        }
    }

    private static class RegionalStackDialogControls {
        final Choice regionSizeChoice;
        final TextField featureSizeField;
        final Checkbox regularizationCheckbox;
        final Checkbox resultsTableCheckbox;

        RegionalStackDialogControls(
                Choice regionSizeChoice,
                TextField featureSizeField,
                Checkbox regularizationCheckbox,
                Checkbox resultsTableCheckbox
        ) {
            this.regionSizeChoice = regionSizeChoice;
            this.featureSizeField = featureSizeField;
            this.regularizationCheckbox = regularizationCheckbox;
            this.resultsTableCheckbox = resultsTableCheckbox;
        }
    }

    private static class SingleImageDialogControls {
        final Choice regionSizeChoice;
        final TextField featureSizeField;
        final Checkbox regularizationCheckbox;
        final Checkbox resultsTableCheckbox;

        SingleImageDialogControls(
                Choice regionSizeChoice,
                TextField featureSizeField,
                Checkbox regularizationCheckbox,
                Checkbox resultsTableCheckbox
        ) {
            this.regionSizeChoice = regionSizeChoice;
            this.featureSizeField = featureSizeField;
            this.regularizationCheckbox = regularizationCheckbox;
            this.resultsTableCheckbox = resultsTableCheckbox;
        }
    }

    private static class RegionData {
        final int number;
        final int gridRow;
        final int gridColumn;
        final int x0;
        final int y0;
        final int width;
        final int height;
        final double centerX;
        final double centerY;

        int bestSlice;
        double maximumFcs;
        boolean reliable;
        double assignedFocusSlice;
        int focusSourceRegion;
        double focusSourceDistance;
        double focusSlice;
        int reconstructionSlice;

        RegionData(
                int number,
                int gridRow,
                int gridColumn,
                int x0,
                int y0,
                int width,
                int height
        ) {
            this.number = number;
            this.gridRow = gridRow;
            this.gridColumn = gridColumn;
            this.x0 = x0;
            this.y0 = y0;
            this.width = width;
            this.height = height;
            this.centerX = x0 + 0.5 * (width - 1);
            this.centerY = y0 + 0.5 * (height - 1);
        }
    }

    private static class Reconstruction {
        final ImageProcessor allInFocus;
        final FloatProcessor focusMap;
        final FloatProcessor qualityMap;
        final ByteProcessor reliabilityMap;

        Reconstruction(
                ImageProcessor allInFocus,
                FloatProcessor focusMap,
                FloatProcessor qualityMap,
                ByteProcessor reliabilityMap
        ) {
            this.allInFocus = allInFocus;
            this.focusMap = focusMap;
            this.qualityMap = qualityMap;
            this.reliabilityMap = reliabilityMap;
        }
    }

    private static class SpectralMeasurement {
        final double crossMagnitude;
        final double powerA;
        final double powerB;

        SpectralMeasurement(double crossMagnitude, double powerA, double powerB) {
            this.crossMagnitude = crossMagnitude;
            this.powerA = powerA;
            this.powerB = powerB;
        }

        double geometricBandPower() {
            if (powerA <= 0.0 || powerB <= 0.0) {
                return 0.0;
            }
            return Math.sqrt(powerA * powerB);
        }
    }

    private static class FftPlan {
        final int splitWidth;
        final int splitHeight;
        final FloatFFT_2D fft;
        final float[] imageA;
        final float[] imageB;
        final int[] includedBins;

        FftPlan(int analysisWidth, int analysisHeight, int featureSize) {
            splitWidth = analysisWidth / 2;
            splitHeight = analysisHeight / 2;
            fft = new FloatFFT_2D(splitHeight, splitWidth);
            imageA = new float[2 * splitWidth * splitHeight];
            imageB = new float[2 * splitWidth * splitHeight];
            includedBins = createHighPassMask(splitWidth, splitHeight, featureSize);
        }

        private static int[] createHighPassMask(
                int splitWidth,
                int splitHeight,
                int featureSize
        ) {
            int[] temporary = new int[splitWidth * splitHeight];
            int count = 0;
            double cutoff = 2.0 / featureSize;
            double cutoffSquared = cutoff * cutoff;

            for (int y = 0; y < splitHeight; y++) {
                int frequencyY = y <= splitHeight / 2 ? y : y - splitHeight;
                double normalizedY = frequencyY / (double) splitHeight;

                for (int x = 0; x < splitWidth; x++) {
                    int frequencyX = x <= splitWidth / 2 ? x : x - splitWidth;
                    double normalizedX = frequencyX / (double) splitWidth;
                    double radiusSquared = normalizedX * normalizedX
                            + normalizedY * normalizedY;

                    if (radiusSquared >= cutoffSquared) {
                        temporary[count++] = y * splitWidth + x;
                    }
                }
            }

            int[] mask = new int[count];
            System.arraycopy(temporary, 0, mask, 0, count);
            return mask;
        }
    }
}
