# Talk Transcript
## COVID-19 Transmission Clusters & Socioeconomic Deprivation in Scotland

*Approximate duration: 7 minutes. Speaking pace: ~140 words per minute.  
Slide cues are marked in bold.*

---

### **[Slide 1 — Title]** *(~30 seconds)*

Thank you for the introduction. I want to talk today about transmission clusters — groups of genomically linked COVID-19 cases — and what they can tell us about how socioeconomic deprivation shaped the spread of the virus across pandemic waves in Scotland. This work uses whole-genome sequencing data to identify clusters and then models what drives how large those clusters get.

---

### **[Slide 2 — Background]** *(~60 seconds)*

So why do clusters matter? When we identify a transmission cluster, we're identifying a group of people who are part of the same chain of transmission — the same outbreak, essentially. A larger cluster means more people infected from a common source, and that's our outcome of interest.

COVID-19 is well-known as a superspreading disease — the majority of transmission is driven by a small minority of events. That has modelling implications I'll come back to.

We used SIMD — the Scottish Index of Multiple Deprivation — as our measure of socioeconomic context. Quintile 1 is the most deprived fifth of Scottish areas, quintile 5 the least deprived. And we had three core questions: does deprivation predict cluster size at all; does the *mixing* of different deprivation groups within a cluster matter; and critically, do these effects change across pandemic waves?

---

### **[Slide 3 — Data & Approach]** *(~60 seconds)*

In terms of the data and method — we worked with 83,889 non-singleton transmission clusters spanning five pandemic waves from Pre-VOC through to Omicron BA.2+.

Two modelling decisions are worth flagging. First, we excluded singleton clusters — those with only one sequenced case. We found their prevalence was very strongly correlated with low sequencing coverage, Spearman r of minus 0.82, which tells us they're mostly a surveillance artefact rather than genuine isolated infections. Second, because sequencing coverage varied enormously over time, we included it as a log-scale offset in the model — that corrects for surveillance intensity without letting it dominate the results.

The model itself is a zero-truncated negative binomial GLM. The confirmed overdispersion parameter of 3.6 tells us superspreading is very much present even after conditioning on all our predictors.

---

### **[Slide 4 — Finding 1: Mixing Dominates]** *(~75 seconds)*

The first and perhaps most striking finding is about socioeconomic mixing *within* clusters. We measured this as the standard deviation of SIMD quintiles across cases in a cluster — so a high value means cases are spread across many different deprivation levels, a low value means they're concentrated in one.

The coefficient on this mixing term was 1.42, with a z-statistic of 129. To put that in practical terms: each unit increase in within-cluster mixing multiplies expected cluster size by about 4.1. You can see that on the right — at zero mixing, Delta-wave Q3 clusters average around 5 sequences; at high cross-quintile mixing, that rises to around 46.

This was by far the strongest predictor in the model — stronger than deprivation level, stronger than epoch. The implication is that it's not just *where* you are in the deprivation distribution that matters, it's whether your social contacts span across it.

---

### **[Slide 5 — Main Result: The Interaction]** *(~90 seconds)*

Now the main result. This figure shows predicted cluster size relative to our reference — Delta wave, SIMD quintile 3 — on the y-axis, with SIMD quintile on the x-axis. Each line is a different pandemic wave, and if the effects were simply additive, these lines would run parallel. They don't.

The most striking feature is Omicron BA.1 in orange. Clusters in the most deprived areas — Q1 — were predicted to be 4.7 times larger than the Delta reference. That's a very large effect. And it's not uniform across quintiles; the interaction is real and statistically significant.

Alpha in red shows essentially the mirror image: less deprived areas saw substantially smaller clusters during Alpha, while deprived areas were relatively spared that protective effect. The least deprived quintile was about 37% smaller than Q1 during Alpha.

Delta and Omicron BA.2+ by contrast show relatively flat lines — a more uniform pattern across deprivation levels.

Pre-VOC is shown for completeness but should be interpreted cautiously — there were only around 2,000 clusters in that epoch and they were geographically concentrated.

---

### **[Slide 6 — Three Distinct Deprivation Gradients]** *(~60 seconds)*

The bar chart makes the three-wave story cleaner to read.

During **Alpha**, there's a clear deprivation gradient — Q1 clusters are the largest, and cluster size falls progressively through to Q5. Deprived communities were disproportionately affected.

During **Delta**, the picture reverses. Q4 and Q5 — less deprived areas — had larger clusters than the middle quintiles. This is consistent with the pattern of behaviour as restrictions lifted in 2021, with indoor social mixing resuming first in more affluent settings.

And during **Omicron BA.1**, everything scales up, but Q1 scales up by far the most. The most deprived communities experienced the largest absolute and relative increase. This represents the clearest health inequality signal across the entire study period.

---

### **[Slide 7 — Conclusions]** *(~45 seconds)*

To summarise. Within-cluster socioeconomic mixing is the dominant driver of cluster size, with a fourfold effect per unit — this suggests that interventions targeting cross-community transmission may be more effective than approaches stratified solely by deprivation level.

The deprivation effect is not fixed — it's wave-specific and requires a full epoch-by-SIMD interaction to capture. Omicron BA.1 produced the starkest inequality we observed, with the most deprived communities experiencing dramatically larger clusters. Alpha showed a steeper deprivation gradient; Delta showed the reverse.

Thank you. Happy to take questions.

---

*Total word count: ~950 words · Estimated duration: 6 min 45 sec at 140 wpm*