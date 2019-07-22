Only mast_50_all_all and mast_50_m30_all gives the correct motif, surprisingly mg_50_30_all doesn't. Also from previous, mg_50_30_1 also gives the correct one. 

Probably redo this with converge, just to see what we'll get. Otherwise report to Igor, using the mg_50_30_1 one, although that's not quite right. But it's k, we can use mast to find the remaining hits. 

Problem is that converge requires seeds, and these are fairly important for converge. What seeds then do we use? Randomly select? Check with Igor on mon, but random selection seems fine. 

The trouble is uniref50's 9k+ are already clusters, so technically we should use that as seeds, and then the 45k one as main. But 9k as seed for converge isn't going to work. 
