Reviewer Commment
=========================================================

Review nuptial gift giving round 2

This is my second round of review on this manuscript. I was enthusiastic during the first round, but had spotted some confusions that I would have liked to see clarified. I have to say that I am a little bit underwhelmed with this new version of the manuscript. While the authors have clearly put efforts into answering some of my numerous comments, I feel like some important concerns on the logic of the model and presentation of the results have been politely left aside. The author’s logic may still be correct (or not), but it is not more explicit in this new version than it was in the first one. I will not give a detailed review here but only insist on those points that I am surprised to see unchanged. This is rather disappointing as I was looking forward to the updated version of an interesting study that I would like to see published.

I did not received the new version without track changes, so I had to open the word document and accept the track changes to work on it. The line numbering I am using in my comment is of the word document with accepted track changes.


Author Response
=========================================================






Reviewer Commment
=========================================================

I have three main concerns. I give an overview here and then detailed comments referring to the
first round of review.

1.) I do not agree with the way that encounter rates are calculated as a function of the sex-ratio in the population. The way it is now, global encounter rates (population density) vary with the sex ratio so that for example when the sex ratio is extremely biased encounter rates become arbitrarily high. My view is that the relative encounter rates when comparing males and females are correct, but the sex-specific encounter rates are wrong. I may be the one who is confused here, but then the logic of how these encounter rates are calculated should be spelled out. If over-estimating encounter rates does not affect the outcome of the model because only the relative encounter rates are important, this is also fine but needs to be spelled out that way. See my point below. I do not think that this is a trivial point, since the evolution of gift searching is likely to modify the time-outs and thus the sex ratio in the mating pool. If encounter rates are not calculated correctly, this could bias the fitness threshold for the evolution of the trait.


Author Response
=========================================================

We greatly appreciate the persistence and the patience of Reviewer 1. We were getting hung up on the interaction rates being consistent with previous literature, but Reviewer 1 is completely correct that only the relative encounter rates were correct while the absolute sex-specific encounter rates are what really matters. Their explanation in their most recent review makes this very clear.

We have revised our model exactly as they have suggested, starting with the absolute encounter rate 'R' instead of 'M' to calculate the absolute rate at which a focal individual encounters a potential mate. Fortunately, the very intuitive thresholds for females and males are preserved; i.e., the theoretical results do not change at all qualitatively. The thresholds are slightly different and linear rather than curvilinear. This is reflected in a new Figure 2. Unfortunately, given this new formulation, there is no closed form solution for beta, but values of beta are easy to calculate through recursion to an arbitrary degree of precision.

Fortunately, no changes to the individual-based model needed to be made because the encounter rates were calculated directly from model output.


Reviewer Commment
=========================================================

2.) I do not understand the discrepancies between the IBM and the analytical model and I feel
like they are being downplayed. For example, the female fitness threshold is qualitatively
different, see my detailed comment and figure below.


Author Response
=========================================================

Again, we are very grateful for Reviewer 1's thoughtful critique. Addressing this concern has helped us reflect on our results and communicate them much more strongly. Discrepancies between the IBM and the analytical model were not intentionally downplayed. We actually now think that they were unintentionally exaggerated by the parameter space that we were illustrating in analytical model (our accidental reversal of threshold colours didn't help). We have resolved this by turning Figure 2 into a two-panel figure, with the second panel now illustrating an area of parameter space that is much more in line with the IBM results of Figure 3.  

Figure 2a now shows the parameter space in the old Figure 2, with the four zones corresponding to four different evolutionary outcomes as before, determined by thresholds for male nuptial gift search (blue) and female choosiness (red). The new Figure 2b shows what happens as conspecific interaction rate (see Reviewer 1's first concern) is decreased from 2 to 1. In this case, zone C (where females would be choosy but males do not search) disappears because the thresholds no longer cross along the range of the x-axis. We are left with three zones: (A) no searching or choosiness, (B) searching but no choosiness, and (D) searching and choosiness. These zones are qualitatively identical to the pattern produced by the IBM in Figure 3, and they are almost exactly predicted by the thresholds shown in Figure 3.

We hope that Reviewer 1 agrees that this new presentation illustrates that the female fitness threshold is not qualitatively different between the analytical model and the IBM. Nevertheless, there are two small quantitative differences for female choosiness that are worth noting. The first is that the female threshold in the IBM (red line of Figure 3) is curved rather than linear. This is because we cannot actually control the interaction rate (R) in the IBM; it isn't fixed as in Figure 2, so we use estimations of the interaction rate to build the threshold line from inequality 5. This threshold line does a good job of predicting when choosiness will evolve, but it is a bit conservative for high nuptial gift search times (alpha; Figure 3). This is likely because the IBM allows for coevolution between male search and female choosiness, which does not exist for the anlaytical model in Figure 2. Note that when female and male traits are fixed while the other evolves (Figure S6.1), the choosiness threshold line appears to be too ambitious -- though this is likely because a much higher sample size is needed for confidence intervals to not overlap zero when parameter values are just over the choosiness threshold in the IBM (simulation time has been a recurring challenge for this manuscript -- each point in Figure 3 represents 3200 replicate simulations, while Figure S6.1 represents 1600).


Reviewer Commment
=========================================================

3.) I think the gamma parameter (increase in offspring relative fitness) needs to reach very high values for gift giving to evolve and I would like to see that discussed (imagine a population genetics or evolutionary dynamics paper where the mutant invading the population has a fitness 6 times higher than the resident). This is very important since the smaller the range of gamma the more differences between the IBM and the analytical model are visible.


Author Response
=========================================================

We agree that Figures 2 and 3 might give readers the misleading impression that very high values of gamma are required for the evolution of male search (and subsequent female choosiness). This is an important point and provides a further rationale for exploring a wider area of parameter space, which the reviewer also suggested. In the supporting information, we now include results from our individual-based model that show an area of parameter space in which male searching evolves at low values of gamma. This complements the analytical model in showing that the evolution of searching is possible given sufficiently low M (now D), even if alpha is low (Figure 2, equation 2). This is even further solidified by our correction of the analytical model brought about by the reviewer's first point, which shows a linear relationship (instead of a curvilinear one) between search time and fitness thresholds -- we're very grateful for these suggestions from the reviewer! We now include two new paragraphs in the Discussion to elaborate on why searching and choosiness are still likely to evolve in some systems, even when gamma is not extremely high (line X). We have also included a whole new section of supporting information which demonstrates that both traits can evolve even when gamma is very low (line X).


Reviewer Commment
=========================================================

4.) New concern: given the discrepancies between the analytical model and IBM as it is, I think that the parameter space should be properly explored. How do Tm,M,Tm, affect the two models? The little arrows on Figure 2 are not enough.


Author Response
=========================================================

We appreciate this new concern, but we hope that our response to Reviewer concerns 1 and 2 demonstrate that there is no major discrepancy between the analytical model and the IBM. Both models indeed illustrate the same qualitative results, so we have provided two separate models leading to identical conclusions about the evolution of nuptial gift searching and choosiness.

Nevertheless, we have taken Reviewer 1's suggestions and tried to explore more of the parameter space for the IBM. A major challenge here has been the computational time needed to run simulations. Even using a compiled programming language (C) and with access to a computing cluster, a full search of parameter space has proved challenging. Inevitably, IBM results will underestimate when gift-giving and choosiness evolve, but evolution should not happen under calculated thresholds.

We have explored how  varying Tm, Tf and M (now R) affected the female and male fitness threshold. Specifically, we re-created figure 3 from the main text under new values of Tm, Tf and M. We explored values both above and below the values used in the main text. 

We found that inequality 2 & 3 of the analytical model do indeed give predictions consistent with calculated thresholds as to changes in Tm, Tf and R. That is, increases in Tm, Tf and R were found to decrease the female fitness threshold and increases in M were found to increase the male fitness threshold, and this is in agreement with the arrows on figure 2 in the main text. We have written a section of supporting information for these results (see Supporting Information S2). With more time, we could run more of these simulations.
























