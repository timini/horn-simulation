# Reading response metrics

Target-band ripple and model feasibility use all original samples and the two requested band edges. They determine whether the candidate satisfies the requested operating band. A narrow notch inside that band cannot be hidden by a broad overall bandwidth.

The secondary `f3` metrics describe the contiguous lobe around the largest sampled peak, at a threshold 3 dB below that peak. Crossings are interpolated between original samples; disconnected lobes are not joined. No extra smoothing or sparse resampling precedes crossing detection.

If the response remains above the threshold at a sweep boundary, an exact cutoff has not been observed. The numeric boundary is retained for compatibility, but JSON also records `f3_low_is_bound` or `f3_high_is_bound`; reports display **≤** for a lower cutoff below the sweep and **≥** for an upper cutoff above the sweep. `bandwidth_is_lower_bound` marks the corresponding lower bound on bandwidth. For example, `f3 high ≥2263 Hz` does not claim that the horn cuts off at 2263 Hz. Extend the sweep to identify that crossing. These bounds are relative to the largest peak within the available sweep; an unseen larger peak can change the threshold.

KPI ripple and mean level describe this observed peak lobe and can differ from the target-band ripple and mean level used for acceptance. Frequency resolution still limits the features that any sampled response can reveal. Physical prediction status is separate from all of these numerical metrics.
