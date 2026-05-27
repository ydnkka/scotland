---
marp: true
theme: roslin
size: 16:9
paginate: true
---
<!-- _class: title -->
<!-- _paginate: false -->

<img src="../logos/uoe/College of Medicine and Veterinary Medicine logo [WHITE].png" alt="College logo" class="title-logo">

<div class="title-layout">

<div class="title-text">

# Identifying Genomic Signatures of SARS-CoV-2 Superspreading Events in Scotland

## Roslin Internal Seminar

<div class="speaker">
Dominic Arthur<br>
PhD student, Kao Group<br>
The Roslin Institute, University of Edinburgh<br>
<b>Supervisors:</b> Prof Rowland R. Kao & Dr Christopher J. Banks
</div>

</div>

<div class="title-figure">
  <img src="images/scottish_map.png" alt="Scottish map">
</div>

</div>

<div class="logo-strip">
  <img src="../logos/uoe/The Roslin Institute logo [WHITE].png" alt="Roslin logo" class="logo">
  <!-- <img src="../logos/uoe/The University of Edinburgh logo [WHITE].png" alt="University of Edinburgh logo" class="logo"> -->
  <img src="../logos/wellcome/wellcome-logo-white.png" alt="Wellcome logo" class="logo">
</div>

<div class="title-date">
28 May 2026
</div>

<!--
Speaker notes:
Good afternoon, everyone. I am Dominic, a final-year PhD student in the Kao cluster. Over the next 20 minutes I want to take you through work that asks whether whether Scotland's COVID-19 genome data can help us detect superspreading events — those situations where one person or one setting drives a disproportionately large fraction of onward transmission These events shaped the COVID-19 pandemic, and the challenge is identifying them from the data we have access to. Most importantly, the genomic diversity of the COVID-19 is quite low for the resolution analysis of interest.
-->

---
<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Research Group

# The Kao group studies how population structure shapes epidemics

<div class="columns">
  <div class="panel">
    <h3>What we study</h3>
    <p>How contact patterns, movement, geography, and host populations change the way infections spread.</p>
    <p>Applications include human pathogens such as SARS-CoV-2 and influenza, and livestock diseases such as foot-and-mouth disease and bovine TB.</p>
  </div>

  <div class="panel">
    <h3>How we study it</h3>
    <p>Contact and movement networks, pathogen genomes, spatial data, and statistical or mechanistic transmission models.</p>
    <p>My project sits in the genomic epidemiology part of that toolkit.</p>
  </div>
</div>

> Collaboration fit: structured contact data, pathogen genomes, surveillance modelling, and spatial epidemiology.

<!--
Speaker notes:
The Kao cluster is interested in how population structure changes disease spread. By population structure I mean things like who comes in contact whom, how animals move between farms, how geographically connected communities are, and how those patterns change over time. The cluster works across both human and animal pathogens — SARS-CoV-2 and influenza on one side, foot-and-mouth disease and bovine TB on the other. Despite the very different biology, the conceptual questions are often the same. My project sits in the genomic epidemiology strand: using virus genome sequences as one source of evidence about how transmission unfolded. If any of those areas sound like something your work connects to — contact data, movement records, pathogen sequences, or surveillance modelling — and are interested in a collaborative project, please come and speak to us afterwards.
-->

---
<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Big question

# What can virus genomes tell us about superspreading?

<div class="columns">
  <div class="panel">
    <h3>What we can infer</h3>
    <p>Closely related genomes sampled close in time point to recent transmission.</p>
    <p>Clusters that are unusually large, grow unusually fast, or seed onward clusters flag where intense transmission may have occurred.</p>
  </div>
  <div class="panel">
    <h3>How we go about it</h3>
    <p>Build clusters from genetic and temporal proximity, then ask which are unusual relative to a null expectation set by sampling effort and epidemic context.</p>
    <p>Treat the standouts as candidates for review.</p>
  </div>
</div>

> Genomes flag where to look. Confirming the setting still needs epidemiological evidence — exposure links, outbreak reports, or contact tracing.

<!--
Speaker notes:
The question I want to put on the table is what genome data can actually tell us about superspreading, and where that inference stops. A superspreading event in the strict epidemiological sense requires a setting and exposure links — something I do not have directly. What I do have is a very large set of virus genomes. Closely related genomes sampled close together in time suggest recent transmission, and some clusters of those genomes are unusually large, grow unusually quickly, or seed onward clusters. Those are the patterns worth flagging.
The harder question is how to do this rigorously: how to define clusters, what null to compare against given uneven sampling and changing epidemic context, and how to be honest about the gap between a genomic candidate and a confirmed event. I will use the word "candidate" deliberately throughout — it means a genomic signal worth reviewing, not a confirmed event.
-->

---
<!-- _class: figure -->
<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

<!-- ## Genomic epidemiology of SARS-CoV-2 in Scotland -->

# 350,000+ genomes, linked to who, where, and when

<div class="figure-panel">
  <img src="images/policy_sequences_over_time.png"
       alt="SARS-CoV-2 sequencing in Scotland">
</div>

<div class="caption">
SARS-CoV-2 genomic surveillance in Scotland, 2020–2023. (A) Daily sequence counts (grey bars) with 7-day smoothed mean (line); coloured bands mark restriction regimes, shaded by intensity (right colour bar). (B) Clade frequency over time (stacked areas, Nextstrain clades); the dotted line shows the proportion of confirmed cases sequenced (right axis). A mean of <strong>19.9%</strong> of cases were sequenced across the period, varying by <strong>15.4 percentage points</strong> between waves.
</div>

<!--
Speaker notes:
This work rests on a rich surveillance resource: PHS and COG-UK generated over 350,000 Scottish SARS-CoV-2 genomes between 2020 and 2023, covering every major variant wave and restriction regime.
What matters as much as the volume is what each genome is linked to. PHS hold individual-level connections between every sequence and the patient's collection date, geography, age, sex, area-level deprivation, and vaccination history. That's what lets us ask not just whether a cluster is unusual, but who is involved and where, against which policy and variant background.
One caveat to flag before the methods: sequencing intensity was not constant — any screen for unusual clusters has to be calibrated against varying sampling effort.
-->

---
<!-- _class: figure -->

<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

# Pipeline clusters similar genomes, connect clusters through time, flag unusual bursts, and model associations

<div class="figure-panel">
  <img src="images/pipeline_schematic.png"
       alt="Pipeline from genomes and metadata to time-windowed genome clusters, linked clusters, candidate flags, and association models">
</div>

> <strong>Clusters</strong> are sets of genetically similar sequences sampled close in time, identified using a method we developed that needs no fixed genetic or duration thresholds.

<!--
Speaker notes:
This is the methods slide I want people to remember. First, I split the epidemic into overlapping time windows and cluster together very similar genomes within each window. Second, I connect related clusters across neighbouring windows, so I can follow how a viral lineage appears to move through the sampled data. Third, I ask which clusters look unusually large or newly expanding compared with other clusters at the same time. The important idea is simple: connect genome clusters through time, then look for unusual bursts.
-->

---
<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Screening logic

# Three questions flag a cluster for review

<div class="columns three">
  <div class="panel">
    <h3>Is it large?</h3>
    <p>Does this cluster contain more sampled infections than its contemporaries in the same time window?</p>
  </div>

  <div class="panel">
    <h3>Is it new?</h3>
    <p>What share of its sequences are newly observed, rather than carried over from an earlier cluster?</p>
  </div>

  <div class="panel">
    <h3>Does it continue?</h3>
    <p>Does the cluster end, carry forward, or seed several later clusters?</p>
  </div>
</div>

> A **candidate** is a cluster that scores unusually on one or more of these questions — a review flag, not a confirmed event.

<!--
Speaker notes:
The screen does not try to prove why any cluster grew. It asks three practical questions, each one framed against contemporaries in the same time window.
First: is it large — does it contain more sampled infections than other clusters from the same period? Second: is it new — what share of its sequences are genuinely newly observed, rather than carried over from an earlier cluster? Third: does it continue — does the signal die out, carry forward, or seed several later clusters?
Clusters need at least six sequences to enter the screen; below that the comparisons are not reliable. A cluster that scores unusually on one or more of these questions becomes a candidate for review. That word matters: candidate, not confirmed event.
-->

---
<!-- _class: figure -->

<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Example signal

# A candidate is a pattern through time, not a named event

<div class="columns">
  <div class="figure-panel">
    <img src="images/fig03_subgraph_AM00027_exemplar.png"
         alt="Example flagged SARS-CoV-2 genome-cluster pattern through time">
  </div>

  <div>
    <h3>How to read it</h3><br>
    <ul>
      <li>Columns are time windows</li><br>
      <li>Circles are genome clusters</li><br>
      <li>Bigger circles contain more sampled genomes</li><br>
      <li>Arrows link genetically related clusters across windows</li><br>
      <li>Coloured circles are flagged for review</li>
    </ul>
  </div>
</div>

<!--
Speaker notes:
This is what one flagged pattern looks like. The coloured circles are the clusters the screen flagged using the three questions from the previous slide — they are unusually large, unusually new, or sit at points where the lineage branches forward into several later clusters.
What this picture is not: a map of a known outbreak setting, or a named event. What it is: the temporal shape of a viral lineage in sampled genome data — enough structure to say this part of the record deserves a closer look.
-->

---
<!-- _class: figure -->

<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Result 1

# Flagged patterns rise and fall across pandemic periods

<div class="figure-panel">
  <img src="images/fig05_candidate_rate_over_time.png"
       alt="Flagged candidate genomic patterns over time by category">
</div>

<div class="caption">A timeline of flagged genomic patterns, not confirmed superspreading events.</div>

<!--
Speaker notes:
When the screen is applied across the surveillance period, the flagged patterns are not evenly distributed through time. Some periods contain more unusual genome clusters than others, and the kind of pattern changes as well. I would be cautious about over-interpreting the peaks. Variant replacement, testing policy, sequencing intensity, and restrictions all affect what appears in the genome data. The useful message is that the screen gives us a structured way to find unusual parts of the epidemic record at national scale.
-->

---
<!-- _class: figure -->

<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Result 2

# What was happening on the Scottish Isles?

<div class="figure-panel">
  <img src="images/fig13_healthboard_map.png"
       alt="Health-board and urban-rural odds ratios for candidate-node membership">
</div>

<div class="caption">
Candidate-node odds ratios by health board (A) and urban–rural class (B); island boards show the strongest enrichment.
</div>

<!--
Speaker notes:
The strongest place signal in the analysis sits in the Scottish islands. Panel A maps the candidate-node odds ratio by health board against Greater Glasgow and Clyde — how much more or less likely sequences from each board are to sit in a flagged cluster. The island boards stand out most: Western Isles at 1.45, Orkney at 1.24, Shetland at 1.14. Panel B shows a more modest urban–rural gradient, with remote rural and remote towns elevated by around 9–11% over large urban areas.
So what was happening there? Honestly, from these data alone, we cannot say. The signal could reflect true transmission dynamics — small, well-connected communities where an introduction spreads rapidly through local networks. But it could equally reflect how surveillance operated in those areas, the timing of introductions relative to local viral diversity, or connectivity factors not captured in the model. The island board populations are also small, which adds real uncertainty to the estimates.
This is exactly the kind of pattern the screen is designed to surface — a geographically specific, falsifiable question to take to other evidence. It is a candidate for external validation, not a conclusion.
-->

---
<!-- _class: figure -->

<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

# Findings need cautious interpretation

<div class="figure-panel">
  <img src="images/validation_bridge.png"
       alt="Bridge from genomic signals to epidemiological validation evidence">
</div>

> Genomes prioritise where to look; epidemiology explains what happened.

<!--
Speaker notes:
This is the conceptual point I want to leave you with. The genomic screen scans years of surveillance data and produces a prioritised review list at national scale — something manual review cannot do. What it cannot do is confirm what happened or why. That requires evidence from outside the graph: setting information, outbreak investigations, contact tracing, and travel or mobility data. Genomics and epidemiology are complementary, not competing. The ideal workflow links them, with genomic flags triggering targeted epidemiological follow-up.
-->

---
<img src="../logos/uoe/The Rosin Institute logo [COLOUR].png" alt="Roslin logo" class="section-logo">

## Conclusion

# A national-scale screen for candidates — and where it leads

<div class="columns">
  <div>
    <h3>What this work shows</h3>
    <ul>
      <li><strong>Detectable signatures</strong> of unusual growth, novelty, and onward spread across 350,000 Scottish genomes.</li><br>
      <li><strong>Structured by place</strong> — candidate membership is most consistently linked to health board, with island boards standing out compared to Glasgow.</li>
    </ul>
  </div>

  <div>
    <h3>Open questions</h3>
    <ul>
      <li>Can setting, mobility, or travel data confirm flagged signals?</li><br>
      <li>How did restrictions, reopening, and testing changes shape what was visible?</li><br>
      <li>Does vaccination alter candidate size, composition, or continuation?</li>
    </ul>
  </div>
</div>

> A scalable, timely screen for candidates — not a replacement for epidemiological investigation.

<!--
Speaker notes:
To bring it together: the contribution is a scalable, retrospective screen for candidate superspreading-like signatures across Scotland's pandemic genome data. Clustering across overlapping time windows and linking clusters into a transition graph identifies nodes with unusual amplification, novelty, and onward spread — at a scale manual review cannot match. The clearest population finding is that candidate membership is structured by place: most consistently by health board, with island boards standing out, and only modestly by urban–rural class.
The open questions are where I'd most welcome input. External validation is the biggest: can setting, mobility, or travel data confirm what the screen flags? The policy layer matters too — restrictions and testing changes shaped which clusters were visible, and separating that from real transmission changes is methodologically important. Vaccination effects on cluster composition are largely unexplored in this framework. And finally, generalisability: does this transfer to influenza or other respiratory pathogens, where genomic surveillance is expanding? If any of these overlap with your work, please come and talk to me.
-->

---
<!-- _class: title -->
<!-- _paginate: false -->

<img src="../logos/uoe/College of Medicine and Veterinary Medicine logo [WHITE].png" alt="College logo" class="title-logo">

<div class="title-layout">

<div class="title-text">

# Acknowledgements

## Thank you

<div class="speaker">
Supervisors: Prof Rowland R. Kao and Dr Christopher J. Banks<br>
Kao Group, The Roslin Institute<br>
Public Health Scotland surveillance and data teams<br>
COG-UK consortium and Scottish sequencing contributors<br>
Wellcome Trust funding
</div>

</div>

<div class="title-figure">
  <img src="images/scottish_map.png" alt="Scottish map">
</div>

</div>

<div class="ack-logo-strip">
  <img src="../logos/uoe/The Roslin Institute logo [WHITE].png" alt="Roslin logo" class="logo">
  <!-- <img src="../logos/uoe/The University of Edinburgh logo [WHITE].png" alt="University of Edinburgh logo" class="logo"> -->
  <img src="../logos/wellcome/wellcome-logo-white.png" alt="Wellcome logo" class="logo">
  <img src="../logos/others/New_COG-UK_logo.jpg" alt="COG-UK logo" class="logo">
  <img src="../logos/others/Public_Health_Scotland_logo.jpg" alt="Public Health Scotland logo" class="logo">
</div>

<!--
Speaker notes:
Thank you very much for listening. I am very grateful to Rowland Kao and Chris Banks for supervision and feedback, the Kao Group for discussion around the network and modelling framing, Public Health Scotland and COG-UK for the surveillance data infrastructure, and Wellcome for funding. I would especially welcome questions or conversations about how to validate the flagged signals, how to interpret the place signal, and whether similar approaches could be useful for other respiratory pathogens.
-->
