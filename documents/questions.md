#### Q: *What should I do if the server crashes?*

**A:** In some ocassional time the web app might disconnect from the server, in that case please reload the page.

<br/>

#### Q: *Why does my file fails to upload into the app?*

**A:** Please make sure that the file format follows the details especified in **Help**. Namely, the data input is a *csv* file with the first column as doses in micromolars​, followed by samples. It is possible to specif header is provided, by default it is true.

<br/>

#### Q: *Can multiple drugs-response pairs be uploaded simmuntaniously?*

**A:** Multiple drugs-responses not necesarilly of the same dimensions are accepted by the app and plotted side by side, please mind saving an identifier for each drug as the first column of the *csv* file.

<br/>

#### Q: *Will the model work on inhibition data (non-decreasing responses) ?*

**A:** Please after uploading the data change the *Inhibition/Viability* switch to "Off"

<br/>

#### Q: *How can I incresae the resolution of my plot?*

**A:** You can change the type of the plots to *.png* or *.jpg* and increase the resolution of the download to a higher number.

<br/>

#### Q: *Is there data available for trying the app?*

**A:** Yes, please feel free to download a example dataset in the following link provided in the **Fit**  &#8594; **Upload** tab.

<br/>

#### Q: *Is there the possibility to change the colors or format of the plot?*

**A:** The colors are picked to be the most distinctived to this models, the source code can be downloaded and manupulited on your own machine. Differente formats for the plots generated please refer to **Options** tab in **ENDS**.

<br/>

#### Q: *How should I cite the data used in the Workshop?*

**A:** For citing the data that was used in the **Workshop** tab inside **Help**  used in the workshop pleas cite this paper ["Intra-tumour diversification in colorectal cancer at the single-cell level"](https://www.nature.com/articles/s41586-018-0024-3).

<br/>

#### Q: *Why is the model fit failing when we select many preprocessing data switches?*

**A:** When outliers are removed and the data is close too being constant the optimizers that solve for the models parameters might fail to converge, or the MH-algorithm that fits the Bayesian model might fail to give good estimates as well. Please feel free to try less preprocessing or download the post-processsed data from **Downloads** and running the codes in your own machine, for example increasing the lenght of the chain in *npB* might produce better estimates. 

<br/>

