

# **Verification and Validation of Acoustic Simulations for Simple Conical Horns**

## **Section 1: The Acoustic Signature of a Conical Horn: A Theoretical Framework**

The verification of any computational model begins with a robust understanding of the underlying physics it aims to represent. For a conical horn loudspeaker, while its geometry appears simple, its acoustic behavior is a product of several interacting physical principles. A correct simulation must accurately capture these phenomena. This section establishes the theoretical foundation for the on-axis and off-axis acoustic performance of a simple conical horn, providing the physical context needed to interpret and validate simulated Sound Pressure Level (SPL) plots.

### **1.1 The Horn as an Acoustic Transformer: Impedance Matching**

A common misconception is that a horn "amplifies" sound. In reality, a horn is a passive acoustic device that does not add energy to the system.1 Its primary function is to act as an acoustic transformer, improving the efficiency of energy transfer between the sound source (the loudspeaker driver) and the surrounding air.2 This principle of impedance matching is fundamental to understanding horn performance.

The diaphragm of a loudspeaker driver is a relatively dense, high-impedance material, while the air into which it radiates is a very light, low-impedance medium. The characteristic acoustic impedance of air is approximately 420 Pa·s/m.1 When a driver attempts to radiate sound directly into this low-impedance air, the mismatch is severe, analogous to trying to propel a boat by pushing oars through the air instead of water.1 This mismatch results in very low efficiency; most of the electrical energy supplied to the driver is converted into heat in the voice coil, with only a small fraction (often less than 1%) being converted into acoustic power.4

The horn mitigates this mismatch. The small cross-sectional area at the "throat" of the horn, where the driver is mounted, restricts the movement of air, presenting a high acoustic impedance to the driver.5 This allows the driver to develop a high sound pressure for a given diaphragm displacement. As the sound wave propagates down the horn's length, the cross-sectional area gradually increases. This gradual flare transforms the acoustic energy from a state of high pressure and low volume velocity at the throat to a state of low pressure and high volume velocity at the "mouth".6 The large area of the mouth is better coupled to the low impedance of the free air, allowing for efficient radiation of sound energy.

The specific acoustic impedance, defined as the ratio of sound pressure (p) to particle velocity (v), or z=p/v, is the key quantity being transformed.1 By creating a more favorable, higher impedance load for the driver to work against, the horn can dramatically increase the overall electroacoustic efficiency, with figures reaching up to 50% in well-designed systems.1 This increased efficiency yields a significant gain in SPL (often 10 dB or more) for the same amplifier power, and it reduces the required diaphragm excursion for a given output level, which in turn lowers mechanical distortion.1

### **1.2 Anatomy of a Conical Horn's On-Axis SPL Response**

The on-axis SPL frequency response of a conical horn is its most recognizable characteristic and a primary target for simulation verification. It is shaped by three distinct regions: the low-frequency roll-off, the passband, and the high-frequency extension.

#### **1.2.1 The Low-Frequency Roll-Off**

A key distinction between horn types is their low-frequency behavior. An exponential horn possesses a well-defined acoustic cutoff frequency (fc​) below which its loading capability drops sharply.5 A simple, straight-walled conical horn, however, does not have a mathematically defined cutoff frequency in the same way.9 Instead, its ability to provide acoustic loading and control the sound radiation is determined by its physical dimensions, primarily the size of its mouth.

The effective low-frequency limit of a conical horn is governed by the relationship between the sound's wavelength and the mouth dimensions. As a widely accepted rule of thumb, a horn begins to radiate efficiently when the circumference of its mouth is approximately equal to one wavelength of the sound frequency.5 Below this frequency, the wavelength becomes too long for the mouth to effectively control, the horn fails to provide adequate impedance loading for the driver, and the acoustic output rolls off. This roll-off is typically more gradual than the sharp cutoff of a textbook exponential horn. Therefore, a valid simulation of a conical horn should not exhibit a steep, filter-like cutoff but rather a progressive decline in output as frequency decreases below this mouth-determined limit.

#### **1.2.2 Passband Characteristics**

The frequency range above the low-frequency roll-off and below the high-frequency limit of the driver is the horn's passband. In an ideal horn, this region would be perfectly flat. In practice, the on-axis SPL of a conical horn is characterized by a series of peaks and nulls, creating a ripple in the frequency response.5

These ripples are primarily the result of reflections within the horn. The most significant reflection occurs at the mouth, where there is an abrupt change in acoustic impedance as the guided wave transitions from the horn walls to free space.7 This impedance discontinuity causes a portion of the wave energy to be reflected back toward the throat. Additional reflections can occur at the throat itself, especially if the transition from the driver to the horn is not perfectly smooth.

The interference between the forward-propagating wave from the driver and these reflected waves creates a pattern of standing waves inside the horn at specific frequencies. This results in acoustic resonances, which manifest as the peaks and dips observed in the SPL plot.10 The frequency, spacing, and magnitude of these ripples are a direct function of the horn's length, throat and mouth dimensions, and flare angle. An accurate simulation must replicate this characteristic ripple pattern, as it is a fundamental acoustic signature of the horn's geometry. Early experimental work by G.W. Stewart in 1920 already documented that the resonance peaks in a conical horn are very closely integral multiples of the fundamental resonance frequency.10

#### **1.2.3 High-Frequency Extension**

The high-frequency performance of a horn loudspeaker system is largely determined by the characteristics of the compression driver itself.9 Most compression drivers exhibit a natural roll-off at high frequencies.12 However, the horn's geometry, particularly at the throat, plays a critical role.

For a simple conical horn that tapers directly from the driver exit, the transition can be abrupt. This can lead to higher-order acoustic modes being generated at the throat, causing irregularities and a premature roll-off in the high-frequency response. Many advanced horn designs incorporate a carefully shaped phase plug or a non-conical throat section (e.g., exponential or hyperbolic) to ensure a smooth transition from the planar or dome-shaped motion of the driver diaphragm to the expanding wavefront in the horn.7 A simulation of a

*simple* conical horn should reflect the potential for these HF issues, and its response will likely differ from a commercial horn that employs such advanced throat geometry.

### **1.3 The Physics of Constant Directivity (and its Limitations)**

Beyond efficiency, the other primary function of a horn is to control the dispersion, or directivity, of the sound.1 Conical horns are the basis for modern "Constant Directivity" (CD) designs.

#### **1.3.1 The "Constant Directivity" (CD) Principle**

The defining feature of a conical horn is its straight-sided walls. These straight walls provide a constant direction of propagation for the sound wave as it expands from the throat to the mouth.6 This geometric simplicity results in a radiation pattern, or coverage angle, that is more consistent across a wide range of frequencies compared to other horn profiles like the exponential horn.14 Exponential horns tend to "beam," meaning their coverage angle narrows significantly as frequency increases, leading to an uneven distribution of sound in the listening area.5 The conical horn's ability to maintain a relatively constant coverage angle is its principal advantage and the reason it forms the core of most modern CD horns.6

#### **1.3.2 The Role of Mouth Size and Frequency**

A horn's ability to control directivity is not infinite; it is fundamentally limited by diffraction. The horn maintains control over its designated coverage angle only for frequencies where the wavelength of sound is smaller than the horn's mouth dimensions.6 As the frequency decreases and the wavelength becomes comparable to or larger than the mouth dimensions, the horn loses its ability to direct the sound, and the waves begin to diffract, or bend, around the edges of the mouth, resulting in a widening coverage pattern.6

This relationship can be approximated by the formula: Θ=K/(d×f), where Θ is the coverage angle, d is the relevant mouth dimension (e.g., width for horizontal coverage), f is the frequency, and K is a constant.16 This equation shows that to maintain a given coverage angle

Θ down to a lower frequency f, the mouth dimension d must be larger. Similarly, to achieve a narrower coverage angle Θ at a given frequency, the mouth dimension d must also be larger. A simulation must correctly model this frequency-dependent loss of pattern control.

#### **1.3.3 The Nuance of "Waistbanding" and Diffraction**

The transition from controlled directivity to wide dispersion is not always simple. Purely conical horns often exhibit a phenomenon known as "waistbanding," where the directivity pattern actually *narrows* for a range of frequencies near the low-frequency limit of its pattern control, before it ultimately widens at lower frequencies.17 This is a measurable artifact that a correct simulation should be able to predict and is a key feature to look for during validation.

Furthermore, the abrupt termination of the horn at the mouth creates complex diffraction effects. While intuition might suggest that diffraction always causes the pattern to widen, the physics are more complex. In some frequency ranges, diffraction can interfere with the direct radiation in a way that causes pattern narrowing.6 To counteract these undesirable effects, many advanced CD horn designs incorporate a secondary, more rapid flare at the mouth. This helps to smooth the impedance transition to free space and reduce diffraction-related anomalies, keeping the directivity more constant.6

It is crucial to recognize that the "simplicity" of a conical horn's geometry is deceptive. Its acoustic behavior is a complex interplay of mouth-size-dependent low-frequency loading, internal reflections causing passband ripple, and directivity characteristics that include both constant-coverage regions and artifacts like waistbanding. Furthermore, most high-performance commercial horns advertised as "Constant Directivity" are not simple conical horns. They are often hybrid designs that couple an exponential or hyperbolic throat section (for superior driver loading) to a conical bell section (for pattern control).7 A simulation of a purely conical horn should, therefore, not be expected to perfectly match the SPL plot of a commercial hybrid horn. This distinction is vital: comparing a simulation to analytical models of a pure conical horn is a direct model-to-model verification, whereas comparing it to a commercial datasheet is a system-level benchmark against a different, more complex geometry. Understanding these nuances is the first step toward a meaningful and accurate validation process.

## **Section 2: Analytical Verification and First-Order Approximations**

Before comparing simulation results to external data, a series of first-order analytical checks can provide a powerful "sanity check." These calculations, based on the fundamental physics of horns, can quickly reveal gross errors in a simulation's setup, geometry, or underlying physical model. If the simulated output deviates significantly from these foundational estimates, it indicates a need to review the model's implementation before proceeding to more detailed validation.

### **2.1 Calculating the Effective Low-Frequency Limit (fc​)**

As established, a conical horn does not have a sharp, well-defined cutoff frequency. However, its performance is fundamentally limited by its mouth size, and it is possible to calculate a practical low-frequency limit, often denoted as fc​, below which the horn ceases to provide efficient loading. A robust check involves triangulating this value using several related rules of thumb.

First, a widely used formula, which is independent of the specific horn flare geometry, relates the required mouth area (SL​) to the desired low-frequency cutoff (fc​).8 This formula is derived from the condition for efficient radiation from a circular piston in an infinite baffle and is given by:

SL​=π1​(2fc​c​)2  
where c is the speed of sound (approximately 344 m/s in air at 20°C). This equation can be rearranged to solve for the estimated cutoff frequency based on the mouth area of the modeled horn:

fc​=2πSL​​c​  
Second, a complementary rule of thumb states that for the horn to effectively guide the wave before it radiates, its axial length (L) should be at least one-quarter of the wavelength (λ) of the lowest frequency of interest.9 This gives another estimate for the cutoff:

fc​≈4Lc​  
Third, a related principle is that the circumference of the mouth (Cmouth​) should be at least one wavelength.5 For a circular mouth of radius

rL​, this yields:

fc​≈2πrL​c​  
These three calculated values for fc​ will not be identical, but they should fall within the same general frequency region. A correctly implemented simulation should show the "knee" of its low-frequency SPL roll-off occurring in the range bounded by these estimates. A significant discrepancy, such as an octave or more, points to a potential error in the model's geometry definition or its boundary conditions at the mouth.

It is critical to avoid applying formulas specific to other horn types. For instance, the cutoff frequency of an infinite exponential horn is given by fc​=m⋅c/(4π), where m is the flare constant.20 While it is possible to approximate a

*local* cutoff frequency for a small segment of a conical horn by treating it as exponential, this formula does not apply to the conical horn as a whole.20 Using it for validation will lead to incorrect conclusions.

### **2.2 The Webster Horn Equation and Its Implications**

The propagation of a one-parameter wave in a rigid-walled horn of varying cross-section, under the assumption of small signal amplitudes, is described by the Webster Horn Equation.22 This is a second-order linear partial differential equation that relates the acoustic pressure

p to the horn's cross-sectional area A(x) along its axis x. For a time-harmonic wave, it can be written as a second-order ordinary differential equation:

dx2d2p​+A(x)1​dxdA(x)​dxdp​+k2p=0  
where k=ω/c is the wave number.

For a conical horn, the area varies as the square of the distance from the virtual apex, A(x)∝x2. The solutions to the Webster equation for this geometry involve spherical Bessel functions, which are not trivial to implement for a quick manual check. However, understanding that this is the governing equation that the simulation pipeline is (or should be) solving is crucial. The simulation is essentially a numerical solver for this equation, coupled with appropriate boundary conditions at the throat and mouth.

More advanced theoretical frameworks recast the Webster equation in terms of forward-propagating (p+​) and backward-propagating (p−​) wave variables.24 This approach provides deeper physical insight into the nature of reflections and energy transfer within the horn. While implementing this wave-variable method is likely beyond the scope of a simple verification check, it underscores the physics that a comprehensive simulation must capture: the forward propagation of the wave from the driver and the generation of a backward-propagating wave due to impedance mismatches, primarily at the mouth.

### **2.3 Throat Impedance: The Driver-Horn Interface**

While the on-axis SPL is the most common output metric, the acoustic impedance at the horn's throat (Zthroat​) is a more fundamental and sensitive property for model validation. Zthroat​ is the acoustic load that the horn presents to the driver diaphragm. It is defined as the ratio of the acoustic pressure to the volume velocity at the throat. Its value is a direct function of the horn's geometry and its radiation into free space, making it an excellent metric for isolating and verifying the horn model itself, separate from the driver model.

A plot of the throat impedance versus frequency is a powerful diagnostic tool. A successful horn design aims to present a high and predominantly resistive (real) acoustic impedance to the driver across its intended passband.8

* **Above the cutoff frequency**, the throat impedance magnitude of a well-behaved horn should be relatively constant and high, and its phase angle should be close to 0 degrees, indicating a resistive load. This condition allows for maximum power transfer from the driver to the air.  
* **Below the cutoff frequency**, the impedance becomes highly reactive, with large positive or negative phase angles, and its magnitude falls dramatically. This indicates that the horn is no longer loading the driver effectively, and most of the energy is reflected rather than radiated.

When comparing horn types, an exponential horn is known to achieve a more constant resistive load closer to its theoretical cutoff frequency. A conical horn's throat impedance typically becomes resistive at a somewhat higher frequency relative to its physical dimensions, indicating that it is generally less efficient at providing low-frequency loading than an exponential horn of similar size.8

If the user's simulation pipeline can calculate and plot the throat impedance, comparing the shape of the magnitude and phase curves to textbook examples provides a rigorous validation of the acoustic model's core calculations. It can reveal subtle errors in the implementation of the wave equation or the radiation boundary conditions at the mouth that might be obscured in the final, system-level SPL plot.

## **Section 3: Empirical Verification: A Curated Database of Reference Designs**

Analytical checks provide a necessary but not sufficient condition for model validation. The next critical step is to compare the simulation's output against reliable empirical data from existing horns. This section provides a curated set of reference data, categorized into DIY projects, commercial products, and published simulations. This database allows for both direct comparison (where geometry is known) and benchmark comparison (where only performance is known).

### **3.1 Public Domain and DIY Project Data**

Do-it-yourself (DIY) and open-source projects offer a unique advantage for validation: the builders often publish complete geometric details alongside their measurement results. This allows for a true end-to-end verification, where the exact reference geometry can be modeled in the custom pipeline and the output compared directly to measured data.

One excellent example is the **Inlow Sound 100 Hz Conical Horn**.25 This is a well-documented project featuring a 12-sided conical horn designed for a 100 Hz low-frequency target, using four 8-inch B\&C 8PE21 drivers. The project website includes a measured on-axis SPL plot, which shows a solid response from approximately 90 Hz to 700 Hz. This serves as a prime reference case because it is a pure conical design with publicly available, real-world measurements.

Another valuable resource is the **Pass Labs "Kleinhorn"** project by Nelson Pass.26 Although this design is an exponential horn, the project documentation is exceptionally thorough. It provides detailed construction plans, driver information, and extensive measurement plots, including frequency response and impedance under various loading and equalization conditions. The methodology and data presentation serve as a best-practice example for how to characterize a horn system.

Online forums, particularly the diyaudio.com community, are a rich source of projects. While the data can be less formally presented, threads often contain detailed discussions, dimensions, and user-generated measurements. For instance, the "Nine-sided conical horn MEH" thread includes a user's comparison of their measurements to simulations, with specific discussion of audible artifacts like midrange narrowing, providing valuable real-world context for these phenomena.17

The variability found in such DIY data should not be seen merely as error but as a valuable signal. An idealized simulation assumes perfect conditions (e.g., rigid walls, perfect geometry, anechoic space). Discrepancies between such a simulation and a real-world DIY measurement are expected and informative. If a DIY measurement shows a stronger low-frequency response than the simulation predicts, it could be due to boundary reinforcement from placing the horn on the floor or in a corner—an effect known as operating in half-space (2π) or quarter-space (π).2 If the measured response is more jagged, it could be due to room reflections not accounted for in the anechoic simulation. By attempting to model these real-world conditions (e.g., by adding a reflective plane to represent a floor), the validation process can evolve into a more advanced model refinement and sensitivity analysis.

To facilitate direct comparison, the following table collates key data from several well-documented public and DIY projects.

**Table 1: DIY and Published Conical Horn Projects \- Dimensions & Performance**

| Project Name/Source | Horn Type | Throat Dimensions | Mouth Dimensions | Length | Driver(s) | Stated fc​ | Link to SPL/Impedance Plot | Notes |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| Inlow Sound Conical Horn | Conical | Singular throat from 4 drivers | Not specified | Not specified | 4x B\&C 8PE21 8" | 100 Hz | 25 | 12-sided construction. Measured SPL shows response from 90 Hz to 700 Hz. |
| Pass Labs Kleinhorn | Exponential (Rear-loaded) | Not specified | \~30 sq ft (\~2.79 m2) | \~9 ft (\~2.74 m) | Lowther DX55 | 30 Hz | 26 | Extensive measurements provided. Serves as a best-practice example for data presentation. |
| DIY Synergy Horn | Conical | 1" throat \+ 4 side ports | Not specified | Not specified | B\&C DE250, 4x Pyle PDMR5 5" | \~230 Hz (mid crossover) | 28 | A complex Synergy design, but demonstrates conical horn principles with measurements. |
| diyAudio "Nine-sided horn" | Conical | 1.4" throat | Not specified | Not specified | Not specified | Not specified | 17 | User notes midrange narrowing, comparing simulation to measurement. |

### **3.2 Commercial Product Benchmarks**

Commercial product datasheets provide a wealth of information for benchmark comparisons. While manufacturers rarely disclose the exact flare profile—many "conical" or "CD" horns are hybrids—the published specifications establish a baseline for what is achievable with professionally engineered systems. By modeling a horn with similar nominal specifications (e.g., coverage angle, throat size), one can validate if the simulation produces results that are in a realistic performance envelope.

Key manufacturers whose datasheets provide useful reference data include:

* **JBL Professional:** Models like the Selenium HL14-50N are explicitly described as conical horns and provide specifications for dispersion and minimum frequency.29 Datasheets for compression drivers like the 2445J often include SPL plots when mounted on specific horns (e.g., the 2380 Bi-Radial horn).31  
* **Beyma:** Offers horns like the TD-164, a 60°x40° constant directivity horn, with a published on-axis frequency response curve that serves as an excellent reference.32 Coaxial drivers like the 10CX300Fe also incorporate a 70° conical horn for the high-frequency section and provide detailed polar and frequency response plots.33  
* **Eighteen Sound:** Provides detailed datasheets for horns like the XT1086 (80°x60° elliptical) and numerous coaxial drivers (e.g., 12NCX750, 10CX650) that feature integrated horns with full performance data.34  
* **Electro-Voice:** A pioneer in constant directivity horns, with extensive documentation. The HR60 and HR4020A datasheets, for example, provide extremely detailed polar plots, beamwidth vs. frequency charts, and on/off-axis response curves.37 These are invaluable for validating directivity simulations.  
* **Community Professional Loudspeakers:** Products like the R.35-3896 and RSH-462 provide detailed specifications including sensitivity, operating range, and nominal beamwidth, setting performance benchmarks for high-output, weather-resistant horn systems.39

The following table summarizes key specifications for a selection of these commercial products, providing a basis for high-level comparison.

**Table 2: Commercial Conical/CD Horns \- Key Specifications & Performance**

| Manufacturer | Model | Type | Throat Dia. | Coverage (H°xV°) | Freq. Range (-10dB) | Sensitivity (dB/1W/1m) | Link to Datasheet/Plot |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| JBL Selenium | HL14-50N | Exponential/Conical | 2.0" | 45° x 45° | Min. Freq: 600 Hz | N/A (Horn only) | 30 |
| Beyma | TD-164 | Constant Directivity | 1.0" | 60° x 40° | Min. Freq: 1200 Hz | N/A (Horn only) | 32 |
| Eighteen Sound | XT1086 | Constant Coverage (Elliptical) | 1.0" | 80° x 60° | Cutoff Freq: 800 Hz | N/A (Horn only) | 35 |
| Electro-Voice | FX12-PRO | Coaxial Conical Horn | 1.4" HF, 12" LF | 80° x 40° | 120 Hz \- 19 kHz | 106 dB (System) | 42 |
| Electro-Voice | HR60 | Constant Directivity | 1.3" | 68° x 41° | Usable from 500 Hz | 103 dB (with driver) | 37 |
| Community | R.35-3896 | Horn-loaded Triaxial | 1" HF, 2x 2.35" MR | 90° x 60° | 80 Hz \- 16 kHz | 98 dB (System) | 39 |

### **3.3 Published Simulation and Research Data**

The most rigorous references come from academic and industry research, where models are often validated against precise experiments.

* **Simulation Software Examples:** Finite Element Method (FEM) and Boundary Element Method (BEM) software packages like COMSOL Multiphysics are the industry standard for high-fidelity acoustic simulation. These companies often provide detailed application notes with full model definitions. COMSOL, for example, has published models for a conical horn lens antenna and a corrugated conical horn, which include complete geometry parameters, material properties, and simulation settings.43 Replicating these specific geometries in a custom pipeline and comparing the results against the output from a trusted commercial solver is a powerful validation technique.  
* **Academic Literature:** Peer-reviewed journals such as the *Journal of the Acoustical Society of America (JASA)* and the *Journal of the Audio Engineering Society (JAES)*, along with AES conference papers, are the ultimate source for validated data. These publications often present analytical, numerical (FEM/BEM), and experimental results for various horn profiles.10 For example, the seminal 1920 paper by G.W. Stewart provides experimental data on the resonance characteristics of conical horns 10, while modern papers present advanced simulation techniques and their experimental validation.46 Accessing these resources can provide high-quality reference cases for direct comparison.

By leveraging this multi-faceted database of DIY, commercial, and academic references, it is possible to build a comprehensive picture of expected conical horn behavior and perform robust comparisons to validate a custom simulation pipeline.

## **Section 4: A Multi-Tiered Workflow for Simulation Validation**

With a firm grasp of the underlying theory and a database of reference designs, a systematic workflow can be employed to validate the custom simulation pipeline. This process is structured in tiers, moving from simple, low-effort checks to more complex and resource-intensive methods. This tiered approach allows for efficient problem diagnosis and builds confidence in the simulation results at each stage.

### **4.1 Tier 1: The Analytical Sanity Check**

This initial tier uses the first-order approximations from Section 2 to perform a quick but essential check for fundamental errors in the model setup.

* **Step 1: Geometry Verification.** The first action is to meticulously verify the implementation of the horn's geometry within the simulation pipeline. This includes the throat area (S1​), mouth area (S2​), and the axial length (L12​). Any error in these primary inputs will invalidate all subsequent results.  
* **Step 2: Low-Frequency fc​ Estimation.** Using the verified mouth area (S2​), calculate the estimated low-frequency limit (fc​) using the analytical formulas presented in Section 2.1. The primary formula is based on mouth area: fc​=c/(2πS2​​).8 Complement this with estimates based on length (  
  fc​≈c/(4L)) and mouth circumference.5  
* **Step 3: Comparison.** Compare these calculated fc​ values with the low-frequency roll-off point on the simulated SPL plot. The \-3 dB or \-6 dB point of the simulated response should fall within the range of the analytical estimates. A significant mismatch (e.g., a factor of two) strongly suggests a fundamental error in the simulation's physics, such as an incorrect radiation impedance model at the mouth or an error in the geometric scaling.

### **4.2 Tier 2: Cross-Simulation with Industry-Standard Tools**

This tier involves comparing the output of the custom pipeline against a widely trusted, publicly available simulation tool. This is arguably the most effective and efficient method for validating the core acoustic model.

* **The Benchmark Tool: Hornresp.** For one-dimensional horn analysis, the freeware program Hornresp, developed by David McBean, is the de facto industry standard for DIY and many professional designers.20 It is a lumped-parameter model that is computationally efficient and has been validated against real-world measurements over many years. Its ability to model conical horns makes it an ideal tool for this task.  
* **Step 1: Model in Hornresp.** Replicate the exact horn geometry in Hornresp. The software has dedicated input fields for a conical horn segment (Con), requiring the throat area (S1), mouth area (S2), and length (L12).20 Ensure the units are consistent.  
* **Step 2: Input Driver Parameters.** Enter the Thiele-Small (T/S) parameters of the driver used in the custom simulation into Hornresp's driver parameter fields. This includes parameters like Fs​, Qts​, Vas​, and Sd​.  
* **Step 3: Simulate and Compare.** Run the Hornresp simulation to generate a predicted on-axis SPL plot. Export this data and overlay it directly with the SPL plot generated by the custom pipeline.  
* **Interpretation.** Since Hornresp and the custom pipeline are attempting to solve the same underlying 1D wave propagation problem, their outputs should be very similar. Minor deviations due to different numerical methods or assumptions are acceptable. However, the overall shape of the SPL curve, the frequency and magnitude of the passband ripples, the location of the low-frequency roll-off, and the overall sensitivity should match closely. A strong correspondence at this stage provides a high degree of confidence that the custom pipeline's acoustic engine is functioning correctly. Numerous tutorials for using Hornresp are available to guide this process.48

### **4.3 Tier 3: Comparative Analysis Against Reference Data**

This tier validates the simulation against the real world by modeling a known, existing horn and comparing the simulated results to its measured performance.

* **Step 1: Select a Reference Case.** From the curated database in Section 3, select a well-documented reference horn. The ideal choice is a pure conical horn with fully specified dimensions and a high-quality measured SPL plot, such as the Inlow Sound 100 Hz horn.25 A commercial horn with a published response curve, like the Beyma TD-164, can also be used, though its exact geometry may need to be estimated.32  
* **Step 2: Model the Reference Case.** Input the dimensions and driver parameters of the chosen *reference horn* into the custom simulation pipeline.  
* **Step 3: Compare Plots.** Overlay the newly generated simulated SPL plot with the actual measured SPL plot from the reference source.  
* **Interpretation.** A perfect one-to-one match is not expected, nor is it the goal. Real-world measurements are influenced by factors not present in an ideal simulation, such as manufacturing tolerances, driver unit variations, and the acoustic environment of the measurement (reflections, boundary loading).17 Instead of seeking a perfect match, the analysis should focus on the correspondence of key features:  
  * Is the simulated low-frequency roll-off at a similar frequency to the measured one?  
  * Are the major resonant peaks and dips in the passband present in both plots, even if their exact amplitudes differ?  
  * Is the overall sensitivity trend across the passband similar?  
  * Is the simulated directivity pattern (if the pipeline is 2D/3D) consistent with the reference horn's specified coverage angle and any published polar data?

The objective is to understand and explain the discrepancies. For example, if the measured bass is higher, it could be due to floor bounce. If the measured high-frequency response is rougher, it could be due to minor imperfections at the horn throat. A successful validation at this tier demonstrates that the simulation can predict the characteristic behavior of a real-world device.

### **4.4 Tier 4 (Optional): The Path to Experimental Validation**

The ultimate form of validation is to build the modeled horn and perform acoustic measurements. This is a resource-intensive process but provides definitive confirmation of the simulation's accuracy.

* **Measurement Hardware and Software.** The modern standard for accessible acoustic measurement involves a calibrated USB microphone (e.g., miniDSP UMIK-1) and powerful analysis software, most notably Room EQ Wizard (REW), which is available for free.51  
* **Measurement Technique.** To validate an anechoic simulation, measurements must be performed in a quasi-anechoic manner. This is achieved by placing the horn and microphone in a large open space (or outdoors) and using an impulse response measurement technique. By applying a time window to the impulse response, it is possible to isolate the direct sound from the horn and discard later-arriving reflections from the ground and other surfaces.51 For low frequencies, where the required time window becomes impractically long, this far-field measurement can be combined with a near-field measurement taken very close to the horn mouth.  
* **Data Processing and Analysis.** The measured data, which should include on-axis and a series of off-axis responses (e.g., in 10° or 15° increments), can be imported into specialized loudspeaker design software like VituixCAD.51 VituixCAD can process this set of polar measurements to generate comprehensive performance plots, including on-axis SPL, listening window average, power response, and directivity index (DI).52 These plots can then be compared directly to the equivalent outputs from the simulation pipeline, providing a complete, multi-faceted validation of the model's predictive capabilities. Detailed guides for this REW-to-VituixCAD workflow are available and widely used in the DIY audio community.53

This tiered workflow provides a scalable and logical path to building confidence in a custom simulation pipeline, starting with simple checks and progressing to full experimental validation as required by the project's goals.

## **Section 5: Synthesis, Interpretation, and Final Recommendations**

The validation process does not end with the generation of comparative plots; it concludes with the interpretation of those results. Discrepancies between simulation and reference data are not necessarily failures of the model but are often valuable sources of information that can lead to a more refined and accurate simulation.

### **5.1 Interpreting Discrepancies: A Diagnostic Guide**

When the simulated SPL plot differs from analytical predictions or reference measurements, the nature of the discrepancy can help diagnose the potential source of error within the custom pipeline.

* **Error in Low-Frequency Cutoff (fc​):** If the simulated low-frequency roll-off point is significantly different from the analytical estimate (Tier 1\) or a reference measurement (Tier 3), the error most likely lies in the model's boundary conditions. Specifically, the acoustic impedance model for the horn mouth radiating into free space may be incorrect. An inaccurate implementation of the horn's mouth area or overall geometric scaling could also be the cause.  
* **Error in Passband Ripple (Frequency or Amplitude):** If the frequencies of the resonant peaks and dips in the passband do not align with reference data, this could point to an error in the simulated length of the horn. If the amplitude of the ripples is incorrect, this suggests an issue with the model's handling of reflections. This could be an incorrect reflection coefficient at the mouth or an inaccurate representation of the acoustic impedance at the throat, which affects how the driver interacts with the reflected wave.  
* **Error in Overall Sensitivity (SPL Level):** If the simulated SPL across the passband is consistently higher or lower than expected, the issue could stem from several sources. An error in the input driver model (e.g., incorrect Thiele-Small parameters) is a common cause. Alternatively, the pipeline's calculation of radiation efficiency or the fundamental impedance transformation from throat to mouth may be flawed.  
* **Error in Directivity / Off-Axis Response:** This is the most complex domain to validate and typically requires a 2D or 3D simulation. If the simulated coverage angle does not match the horn's design target or reference polar plots (e.g., from an Electro-Voice datasheet 37), the error is likely in the model's handling of wave propagation across the wavefront or its calculation of diffraction at the mouth. One-dimensional models like Hornresp cannot predict off-axis performance, so this level of validation requires a more sophisticated simulation engine. Comparing simulated polar plots to those generated from real measurements in a tool like VituixCAD provides the most rigorous test.54

### **5.2 A Unified Strategy for Validation**

To effectively verify that the plots from a custom-built pipeline are "vaguely correct," a unified, multi-tiered strategy is recommended. This approach balances effort with the level of confidence achieved at each stage.

1. **Foundation (Tier 1):** Begin with the non-negotiable analytical sanity checks. Verify the model's geometry and calculate the expected low-frequency limit (fc​) using the formulas based on mouth area and horn length. This is a low-effort, high-reward step that catches fundamental errors immediately.  
2. **Core Validation (Tier 2):** The most powerful and efficient step for validating the core acoustic engine is to perform a cross-simulation using Hornresp. Given that Hornresp is a widely trusted 1D solver, achieving a close match between its output and the custom pipeline's output for the same horn and driver provides a very high degree of confidence in the fundamental physics implementation. For the stated goal of the query, a successful Tier 2 validation may be sufficient.  
3. **Benchmarking (Tier 3):** Use the curated reference data of DIY and commercial horns to benchmark the simulation's performance against the real world. This step provides crucial context. It helps in understanding how an idealized simple conical horn behaves relative to more complex hybrid designs and how real-world factors like room boundaries and construction methods affect the measured response.  
4. **Ultimate Proof (Tier 4):** If the goal is academic publication, commercial product development, or absolute certainty in the model's predictive power, then experimental validation is necessary. This involves building the horn and conducting quasi-anechoic measurements, processing the data with tools like REW and VituixCAD, and performing a direct, feature-by-feature comparison with the simulation's output.

In conclusion, for the purpose of verifying that a custom simulation's SPL plots are broadly accurate, the recommended path is to proceed sequentially through Tiers 1 and 2\. A successful analytical check followed by a strong correlation with the results from Hornresp would provide robust evidence that the custom pipeline is correctly modeling the primary acoustic behavior of a simple conical horn. The reference data in Tier 3 should then be used to understand how the idealized simulation relates to the performance of practical, real-world loudspeaker systems.

#### **Works cited**

1. Acoustical horns and id waveguides, accessed on July 7, 2025, [https://www.rintelen.ch/download/JMMLC\_horns\_lecture\_etf10.pdf](https://www.rintelen.ch/download/JMMLC_horns_lecture_etf10.pdf)  
2. Horn (acoustic) \- Wikipedia, accessed on July 7, 2025, [https://en.wikipedia.org/wiki/Horn\_(acoustic)](https://en.wikipedia.org/wiki/Horn_\(acoustic\))  
3. 7\. The World Through Sound: Acoustic Impedance, accessed on July 7, 2025, [https://acousticstoday.org/7-the-world-through-sound-acoustic-impedance/](https://acousticstoday.org/7-the-world-through-sound-acoustic-impedance/)  
4. Horn Design \- J H S Audio, accessed on July 7, 2025, [https://www.jhsaudio.com/design.html](https://www.jhsaudio.com/design.html)  
5. Horn loudspeaker \- Wikipedia, accessed on July 7, 2025, [https://en.wikipedia.org/wiki/Horn\_loudspeaker](https://en.wikipedia.org/wiki/Horn_loudspeaker)  
6. Understanding Horns, Part 2 | FOH | Front of House Magazine, accessed on July 7, 2025, [https://fohonline.com/articles/speaking-of-speakers/understanding-horns-part-2/](https://fohonline.com/articles/speaking-of-speakers/understanding-horns-part-2/)  
7. Loudspeaker Horns \- A Crash Course | FOH | Front of House Magazine, accessed on July 7, 2025, [https://fohonline.com/articles/tech-feature/loudspeaker-horns-a-crash-course/](https://fohonline.com/articles/tech-feature/loudspeaker-horns-a-crash-course/)  
8. Section 5.0 : Horn Physics \- Quarter Wavelength Loudspeaker Design, accessed on July 7, 2025, [http://www.quarter-wave.com/Horns/Horn\_Physics.pdf](http://www.quarter-wave.com/Horns/Horn_Physics.pdf)  
9. A few of many questions about horn design. \- The Klipsch Audio Community, accessed on July 7, 2025, [https://community.klipsch.com/index.php?/topic/133770-a-few-of-many-questions-about-horn-design/](https://community.klipsch.com/index.php?/topic/133770-a-few-of-many-questions-about-horn-design/)  
10. The Performance of Conical Horns | Phys. Rev., accessed on July 7, 2025, [https://link.aps.org/doi/10.1103/PhysRev.16.313](https://link.aps.org/doi/10.1103/PhysRev.16.313)  
11. THE PERFORMANCE OF CONICAL HORNS. \- Zenodo, accessed on July 7, 2025, [https://zenodo.org/record/2272629/files/article.pdf](https://zenodo.org/record/2272629/files/article.pdf)  
12. Horn Speakers \- Is it me or....... \- Audio Science Review (ASR) Forum, accessed on July 7, 2025, [https://www.audiosciencereview.com/forum/index.php?threads/horn-speakers-is-it-me-or.9633/page-25](https://www.audiosciencereview.com/forum/index.php?threads/horn-speakers-is-it-me-or.9633/page-25)  
13. Influence of driver shape (dome, cone) on sound | Audio Science Review (ASR) Forum, accessed on July 7, 2025, [https://www.audiosciencereview.com/forum/index.php?threads/influence-of-driver-shape-dome-cone-on-sound.9623/](https://www.audiosciencereview.com/forum/index.php?threads/influence-of-driver-shape-dome-cone-on-sound.9623/)  
14. Constant Directivity Horns: A horn provides more sound pressure level (SPL) at a given listening area by increasing the directiv, accessed on July 7, 2025, [https://timan.cs.illinois.edu/pdfjs\_csintro//static/slides/PHYS406/P406POM\_Lect10/PHYS406---P406POM\_Lect10---slide14.pdf](https://timan.cs.illinois.edu/pdfjs_csintro//static/slides/PHYS406/P406POM_Lect10/PHYS406---P406POM_Lect10---slide14.pdf)  
15. Conical horn loudspeaker design : r/diyaudio \- Reddit, accessed on July 7, 2025, [https://www.reddit.com/r/diyaudio/comments/8yzm59/conical\_horn\_loudspeaker\_design/](https://www.reddit.com/r/diyaudio/comments/8yzm59/conical_horn_loudspeaker_design/)  
16. Understanding Horn Directivity Control \- audioXpress, accessed on July 7, 2025, [https://audioxpress.com/article/Understanding-Horn-Directivity-Control](https://audioxpress.com/article/Understanding-Horn-Directivity-Control)  
17. Nine-sided conical horn MEH | Page 6 \- diyAudio, accessed on July 7, 2025, [https://www.diyaudio.com/community/threads/nine-sided-conical-horn-meh.407394/page-6](https://www.diyaudio.com/community/threads/nine-sided-conical-horn-meh.407394/page-6)  
18. Synergy Calc V5, accessed on July 7, 2025, [http://libinst.com/SynergyCalc/Synergy%20Calc%20V5.pdf](http://libinst.com/SynergyCalc/Synergy%20Calc%20V5.pdf)  
19. The Quadratic-Throat Waveguide®: | Peavey Commercial Audio, accessed on July 7, 2025, [https://peaveycommercialaudio.com/wp-content/uploads/2020/06/Quadratic-Throat-Waveguide-by-Charles-E-Hughes.pdf](https://peaveycommercialaudio.com/wp-content/uploads/2020/06/Quadratic-Throat-Waveguide-by-Charles-E-Hughes.pdf)  
20. Hornresp | Page 589 | diyAudio, accessed on July 7, 2025, [https://www.diyaudio.com/community/threads/hornresp.119854/post-6546716](https://www.diyaudio.com/community/threads/hornresp.119854/post-6546716)  
21. Fc of a horn \- The Klipsch Audio Community, accessed on July 7, 2025, [https://community.klipsch.com/topic/238875-fc-of-a-horn/](https://community.klipsch.com/topic/238875-fc-of-a-horn/)  
22. The Sturm-Louville-Webster Horn equation \- Index of /, accessed on July 7, 2025, [https://jontalle.web.engr.illinois.edu/uploads/403/Horns.pdf](https://jontalle.web.engr.illinois.edu/uploads/403/Horns.pdf)  
23. Horn Theory: An Introduction, Part 1 \- audioXpress, accessed on July 7, 2025, [https://audioxpress.com/assets/upload/files/kolbrek2884.pdf](https://audioxpress.com/assets/upload/files/kolbrek2884.pdf)  
24. Reflectance of acoustic horns and solution of the inverse problem \- PMC \- PubMed Central, accessed on July 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3316681/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3316681/)  
25. 100 hz 12 sided conical horn \- The Paper Horn by Inlow Sound, accessed on July 7, 2025, [https://inlowsound.weebly.com/100-hz-12-sided-conical-horn.html](https://inlowsound.weebly.com/100-hz-12-sided-conical-horn.html)  
26. The Kleinhorn Part 1 \- PassDiy, accessed on July 7, 2025, [https://www.passdiy.com/project/speakers/the-kleinhorn-part-1](https://www.passdiy.com/project/speakers/the-kleinhorn-part-1)  
27. The Kleinhorn Part 2 \- PassDiy, accessed on July 7, 2025, [https://www.passdiy.com/project/speakers/the-kleinhorn-part-2](https://www.passdiy.com/project/speakers/the-kleinhorn-part-2)  
28. DIY point source horn \- Red Spade Audio, accessed on July 7, 2025, [http://redspade-audio.blogspot.com/2011/03/diy-synergy-horn.html](http://redspade-audio.blogspot.com/2011/03/diy-synergy-horn.html)  
29. Order JBL Selenium HL14-50N Exponential Horn \- SoundImports, accessed on July 7, 2025, [https://www.soundimports.eu/en/jbl-selenium-hl14-50n.html](https://www.soundimports.eu/en/jbl-selenium-hl14-50n.html)  
30. JBL Selenium HL14-50N 2" Exponential Horn 45x45 4-Bolt \- Parts Express, accessed on July 7, 2025, [https://www.parts-express.com/Selenium-HL14-50N-2-Exponential-Horn-45x45-4-Bolt-264-316](https://www.parts-express.com/Selenium-HL14-50N-2-Exponential-Horn-45x45-4-Bolt-264-316)  
31. 2445J Information \- JBL Professional, accessed on July 7, 2025, [https://jblpro.com/en-US/site\_elements/2445j-information](https://jblpro.com/en-US/site_elements/2445j-information)  
32. TD-164, accessed on July 7, 2025, [https://www.beyma.com/speakers/Fichas\_Tecnicas/beyma-speakers-data-sheet-horn-TD164.pdf](https://www.beyma.com/speakers/Fichas_Tecnicas/beyma-speakers-data-sheet-horn-TD164.pdf)  
33. Coaxial speaker Beyma 10CX300Fe, 8+16 ohm, 10 inch, accessed on July 7, 2025, [https://en.toutlehautparleur.com/coaxial-speaker-beyma-10cx300fe-8-16-ohm-10-inch-3.html](https://en.toutlehautparleur.com/coaxial-speaker-beyma-10cx300fe-8-16-ohm-10-inch-3.html)  
34. 12NCX750 \- Eighteen Sound \- Professional loudspeakers, accessed on July 7, 2025, [https://www.eighteensound.it/en/products/coaxial/12-0/8/12NCX750](https://www.eighteensound.it/en/products/coaxial/12-0/8/12NCX750)  
35. XT1086 \- Eighteen Sound \- Professional loudspeakers, accessed on July 7, 2025, [https://www.eighteensound.it/en/products/horn/1-0/0/XT1086](https://www.eighteensound.it/en/products/horn/1-0/0/XT1086)  
36. 10CX650 \- Eighteen Sound \- Professional loudspeakers, accessed on July 7, 2025, [https://eighteensound.it/en/products/coaxial/10-0/8/10CX650](https://eighteensound.it/en/products/coaxial/10-0/8/10CX650)  
37. HR60 EDS.pdf \- Electro-Voice, accessed on July 7, 2025, [https://products.electrovoice.com/binary/HR60%20EDS.pdf](https://products.electrovoice.com/binary/HR60%20EDS.pdf)  
38. Electro-Voice HR4020A, accessed on July 7, 2025, [https://products.electrovoice.com/binary/HR4020A%20EDS.pdf](https://products.electrovoice.com/binary/HR4020A%20EDS.pdf)  
39. DATA SHEET \- COMMUNITY R Series \- R.35-3896 THREE-WAY FULL-RANGE (90° X 60°) \- Biamp, accessed on July 7, 2025, [https://downloads.biamp.com/assets/docs/default-source/data-sheets/biamp\_data\_sheet\_community\_r-35-3896\_oct21.pdf?sfvrsn=9d51cc\_8](https://downloads.biamp.com/assets/docs/default-source/data-sheets/biamp_data_sheet_community_r-35-3896_oct21.pdf?sfvrsn=9d51cc_8)  
40. Community RSH-462 Horn Speaker \- Sound Productions, accessed on July 7, 2025, [https://www.soundpro.com/community-rsh-462-horn-speaker/](https://www.soundpro.com/community-rsh-462-horn-speaker/)  
41. Community R.35-3896 8" 3-Way Horn-Loaded Weather Resistant Loudspeaker, accessed on July 7, 2025, [https://www.proacousticsusa.com/community-r-35-3896-8-inch-three-way-horn-loaded-loudspeaker.html](https://www.proacousticsusa.com/community-r-35-3896-8-inch-three-way-horn-loaded-loudspeaker.html)  
42. FX12-PRO \- Products | Dynacord, accessed on July 7, 2025, [https://products.dynacord.com/downloadfile.php?id=975026](https://products.dynacord.com/downloadfile.php?id=975026)  
43. Fast Numerical Modeling of a Conical Horn Lens Antenna \- COMSOL, accessed on July 7, 2025, [https://www.comsol.com/model/download/950551/models.rf.conical\_horn\_lens\_antenna.pdf](https://www.comsol.com/model/download/950551/models.rf.conical_horn_lens_antenna.pdf)  
44. Corrugated Circular Horn Antenna Simulator \- COMSOL Multiphysics Application Library, accessed on July 7, 2025, [https://www.comsol.com/model/download/1236311/applications.corrugated\_circular\_horn\_antenna.pdf](https://www.comsol.com/model/download/1236311/applications.corrugated_circular_horn_antenna.pdf)  
45. Reduction of distortion in conical horn loudspeakers at high levels \- ResearchGate, accessed on July 7, 2025, [https://www.researchgate.net/publication/267551388\_Reduction\_of\_distortion\_in\_conical\_horn\_loudspeakers\_at\_high\_levels](https://www.researchgate.net/publication/267551388_Reduction_of_distortion_in_conical_horn_loudspeakers_at_high_levels)  
46. A MODELING AND MEASUREMENT STUDY OF ACOUSTIC HORNS by John T. Post and Elmer L. Hixson May 1994 Electroacoustics Research Labora \- AudioRoundTable.com, accessed on July 7, 2025, [https://audioroundtable.com/misc/post\_hixson\_horns.pdf](https://audioroundtable.com/misc/post_hixson_horns.pdf)  
47. One-Dimensional Acoustic Models of Horns and Comparison with Measurements \- ResearchGate, accessed on July 7, 2025, [https://www.researchgate.net/profile/Thomas-Helie/publication/262853279\_One-Dimensional\_Acoustic\_Models\_of\_Horns\_and\_Comparison\_with\_Measurements/links/5fb13b3a299bf10c36809b04/One-Dimensional-Acoustic-Models-of-Horns-and-Comparison-with-Measurements.pdf](https://www.researchgate.net/profile/Thomas-Helie/publication/262853279_One-Dimensional_Acoustic_Models_of_Horns_and_Comparison_with_Measurements/links/5fb13b3a299bf10c36809b04/One-Dimensional-Acoustic-Models-of-Horns-and-Comparison-with-Measurements.pdf)  
48. Horn Speaker Design Tutorial \- Start Here If You're New \- YouTube, accessed on July 7, 2025, [https://www.youtube.com/watch?v=sENqoDWBq\_c](https://www.youtube.com/watch?v=sENqoDWBq_c)  
49. HORNRESP TUTORIAL Don Radick Version Date Notes 1.0 2-3-25 Initial release 2.0 2-6-25 Added tapered TL example., accessed on July 7, 2025, [https://diy.midwestaudio.club/uploads/editor/qg/s3id8gnkucwe.pdf](https://diy.midwestaudio.club/uploads/editor/qg/s3id8gnkucwe.pdf)  
50. Hornresp tutorial Part 3 \- offset design, a bass reflex and introduction to more terminology, accessed on July 7, 2025, [https://www.youtube.com/watch?v=e3EuVpJnaVY](https://www.youtube.com/watch?v=e3EuVpJnaVY)  
51. VituixCAD help 2.0, accessed on July 7, 2025, [https://kimmosaunisto.net/Software/VituixCAD/VituixCAD\_help\_20.html](https://kimmosaunisto.net/Software/VituixCAD/VituixCAD_help_20.html)  
52. Isobaric And Polar Plotting With miniDSP, REW And VituixCAD, accessed on July 7, 2025, [https://www.minidsp.com/applications/rew/415-isobaric-and-polar-plotting-with-minidsp-rew-and-vituixcad](https://www.minidsp.com/applications/rew/415-isobaric-and-polar-plotting-with-minidsp-rew-and-vituixcad)  
53. Using A Contour Plot for a Polar Plot for Directivity Matching | Audio Science Review (ASR) Forum, accessed on July 7, 2025, [https://www.audiosciencereview.com/forum/index.php?threads/using-a-contour-plot-for-a-polar-plot-for-directivity-matching.45596/](https://www.audiosciencereview.com/forum/index.php?threads/using-a-contour-plot-for-a-polar-plot-for-directivity-matching.45596/)  
54. Vertical Directivity, Grabbing the bull by the horns \- Vandermill-Audio, accessed on July 7, 2025, [https://www.vandermill-audio.nl/vertical-directivity-grabbing-the-bull-by-the-horns/](https://www.vandermill-audio.nl/vertical-directivity-grabbing-the-bull-by-the-horns/)  
55. Measurements with Room EQ Wizard (REW) for crossover simulation with VituixCAD 2, accessed on July 7, 2025, [https://kimmosaunisto.net/Software/VituixCAD/VituixCAD\_Measurement\_REW.pdf](https://kimmosaunisto.net/Software/VituixCAD/VituixCAD_Measurement_REW.pdf)