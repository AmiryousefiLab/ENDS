#### Q: *What should I do if the server crashes?*

**A:** Due to internal Shiny Servers structures, sometimes the web app might momentarily disconnect from the server. In these cases please reload the page and restart the analysis.

<br/>

#### Q: *Why does my file fails to upload into the app?*

**A:** Please make sure that the file format follows the details specified in **Help**. Namely, the data input is a *csv* file with the first column as doses in micromolars​, followed by samples. It is possible to specific header is provided, by default it is true.

<br/>

#### Q: *How can I increase the resolution of my plot?*

**A:** You can change the type of the plots to *.png* or *.jpg* and increase the resolution of the download to a higher number.

<br/>

#### Q: *Is there a data set available for trying the app?*

**A:** Yes, please feel free to download an example dataset in the following link provided in the **Fit**  &#8594; **Upload** tab.

<br/>

#### Q: *Is there the possibility to change the colours or format of the plot?*

**A:** The colours are picked to be the most distinctive for each models, the source code can be downloaded and manipulated on your own machine. For different formats for the plots generated please refer to **Options** tab in **ENDS**.

<br/>

#### Q: *How should I cite the data used in the Workshop?*

**A:** For citing the data that was used in the **Workshop** tab inside **Help**  used in the workshop pleas cite this paper ["Intra-tumour diversification in colorectal cancer at the single-cell level"](https://www.nature.com/articles/s41586-018-0024-3).

<br/>

#### Q: *Why is the model fit failing when we select many preprocessing data switches?*

**A:** When outliers are removed and the data is close too being constant the optimizers that solve for the models parameters might fail to converge, or the MH-algorithm that fits the Bayesian model might fail to give good estimates as well. Please feel free to try less preprocessing or download the post-processed data from **Downloads** and running the codes in your own machine, for example increasing the length of the chain in *npB* might produce better estimates. 

<br/>

