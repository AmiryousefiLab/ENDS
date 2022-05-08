#### Q: *What should I do if the server crashes?*

**A:** Occasionally, the web app might disconnect from the server, in that case please reload the page.

<br/>

#### Q: *Why does my file fails to upload into the app?*

**A:** Please make sure that the file format follows the details specified in **Help**. Namely, the data input is a *csv* file with the first column as doses in micromolars​, followed by samples. It is possible to specify the header and that is by default TRUE.

<br/>

#### Q: *Can multiple drugs-response pairs be uploaded simultaneously?*

**A:** Multiple drugs-responses not necessarily of the same dimensions are accepted by the app and plotted side by side, please mind saving an identifier for each drug as the first column of the *csv* file.

<br/>

#### Q: *Will the model work on inhibition data (non-decreasing responses) ?*

**A:** The web app accepts both viability and inhibition data, after uploading the data select *Inhibition/Viability* switch to "Off" for Inhibition and "On" for Viability. 

<br/>

#### Q: *How can I increase the resolution of my plot?*

**A:** You can change the type of the plots to *.png* or *.jpg* and increase the resolution of the download to a higher number.

<br/>

#### Q: *The x-axis is too jammed and is not easy to read?*

**A:** Probably you have numerous concentrations in the drug-dose range. Please turn off the Concentration option in the plotting options. 

<br/>

#### Q: *Is there any data available for trying the app?*

**A:** Yes, please feel free to download an example dataset in the following link provided in the **Fit**  &#8594; **Upload** tab.

<br/>

#### Q: *Is there the possibility to change the colours or format of the plot?*

**A:** The colours are picked to be the most distinctive to this models, the source code can be downloaded and manipulated on your own machine. Different formats for the plots generated please refer to **Options** tab in **ENDS**.

<br/>

#### Q: *How should I cite the data used in the Workshop?*

**A:** For citing the data that was used in the **Workshop** tab inside **Help**  used in the workshop pleas cite this paper ["Intra-tumour diversification in colorectal cancer at the single-cell level"](https://www.nature.com/articles/s41586-018-0024-3).

<br/>

#### Q: *Why is the model fit failing when we select many preprocessing data switches?*

**A:** When outliers are removed and the data is close too being constant the optimisers that solve for the models parameters might fail to converge, or the MH-algorithm that fits the Bayesian model might fail to give good estimates as well. Please feel free to try less preprocessing or download the post-processed data from **Downloads** and running the codes in your own machine, for example increasing the length of the chain in *npB* might produce better estimates. 

<br/>

#### Q: *Epistemic?*

**A:** Yes. Epistemic, as derived from the Greek word *“επιστήμη”* (as often opposed to aleatory uncertainty), is colloquially defined as “…knowledge or relating to knowledge”. However, on a more formal level, it is related to the sureness to the maximum factual level possible such that all the left error or uncertainty is reflected on the data available about a phenomenon. On the other hand, the aleatory type of uncertainty stems not only from the data but also embrace the uncertainty incurred due to our judgment of a phenomenon, so it has more of the subjectivity flavour into that. For example if one counts the number of 6 red balls out of 10 in a bag can say “the probability of a red ball in the bag is 0.6”, while a person who have got the sample of 2 and sees one of them being red get the “the probability of red is 0.5" which is not necessarily the fact but also have the modelling (generalisation) error in its formation and hence is of lower merit compared to the former case. Also note that, not quite completely, but these two concepts often have been paralleled with the objective and subjective judgement. Our emphasis on the nonparametric models here affirm the fact about our position that, in the drug-response settings (specially with low number of repetitions, systematic sampling of the x-axis, etc.), one should be as objective as possible as the proposed parametric models often have a redundant level of aleatory uncertainty that is making the downstream inferential task non-reliable, while on the other hand the simpler (epistemic) nonparametric models are bereft of these pitfalls.

<br/>
