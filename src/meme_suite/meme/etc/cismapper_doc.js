//
// cismapper-doc.js
//
// Function to replace the innerHTML of element "id" with the HTML indicated by "doc_type".
// If "id" is the empty string, the HTML text is just returned.
//
function print_cismapper_doc(id, doc_type) {
  var html;
  switch (doc_type) {
    case 'html-file-short':
      html = `
	CisMapper outputs an HTML file that provides the results in a human-readable format; 
	this file allows interactive selection, filtering and sorting of the potential regulatory links. 
      `;
      break;
    case 'html-file':
      html = `<p>` + print_cismapper_doc("", 'html-file-short') + `</p>`;
      break;
    case 'tsv-general':
      html = `
	The first line in the file contains the (tab-separated) names of the fields.
	Your command line and other program information is given at the end of the file in 
	comment lines starting with the character '#'.
	The names and meanings of each of the fields are described in the table below.
      `;
      break;
    case 'gene-targets-tsv-short':
      html = `
	CisMapper outputs a gene-centric tab-separated values (TSV) file ('gene_targets.tsv') 
	that contains one line for each gene for which a potential regulatory link was found.
	The lines are sorted in order of increasing Gene_Score.
      `;
      break;
    case 'gene-targets-tsv':
      html = `
	<p>` +
          print_cismapper_doc("", 'gene-targets-tsv-short') + 
	  print_cismapper_doc("", 'tsv-general') + 
          get_cismapper_doc_text('gene-table') + `
	</p>`;
      break;
    case 'gene-elements-tsv-short':
      html = `
	CisMapper outputs an RE-centric tab-separated values (TSV) file ('gene_elements.tsv') 
	that contains one line for each potential regulatory link that was found.
	The lines are grouped by Gene, and the groups are sorted in order of increasing Gene_Score.
      `;
      break;
    case 'gene-elements-tsv':
      html = `
	<p>` + 
          print_cismapper_doc("", 'gene-elements-tsv-short') + 
	  print_cismapper_doc("", 'tsv-general') + 
          get_cismapper_doc_text('gene-table') + `
	</p>`;
      break;
    case 'tss-targets-tsv-short':
      html = `
	If your annotation file provides TSS information,
	CisMapper outputs a TSS-centric tab-separated values (TSV) file ('tss_targets.tsv') 
	that contains one line for each TSS for which a potential regulatory link was found.
	The lines are sorted in order of increasing TSS_Score.
      `;
      break;
    case 'tss-targets-tsv':
      html = `
	<p>` + 
          print_cismapper_doc("", 'tss-targets-tsv-short') + 
	  print_cismapper_doc("", 'tsv-general') + 
          get_cismapper_doc_text('tss-table') + `
	</p>`;
      break;
    case 'tss-elements-tsv-short':
      html = `
          If your annotation file provides TSS information,
	  CisMapper outputs an RE-centric tab-separated values (TSV) file ('tss_elements.tsv') 
	  that contains one line for each potential regulatory link that was found.
	  The lines are grouped by TSS, and the groups are sorted in order of increasing TSS_Score.
	  Within groups, the lines are sorted first by Score, then by TSS, and then by RE_Locus.
      `;
      break;
    case 'tss-elements-tsv':
      html = `
	<p>` + 
          print_cismapper_doc("", 'tss-elements-tsv-short') + 
	  print_cismapper_doc("", 'tsv-general') + 
          get_cismapper_doc_text('tss-table') + `
	</p>`;
      break;
    case 'links-tsv-short':
      html = `
	CisMapper outputs a comprehensive tab-separated values (TSV) file ('links.tsv') 
	that contains one line for each gene for which a potential regulatory link was found.
	The lines are sorted in order of increasing Gene_Score.
      `;
      break;
    case 'links-tsv':
      html = `
	<p>` + 
          print_cismapper_doc("", 'links-tsv-short') + 
	  print_cismapper_doc("", 'tsv-general') + 
	  get_cismapper_doc_text('links-table') + `
	</p>`;
      break;
    case 'tracks-bed-short':
      html = `
	CisMapper outputs a BED file ('tracks.bed') suitable for uploading to the UCSC Genome
	browser for visualizing the potential regulatory links.
      `;
      break;
    case 'tracks-bed':
      html = `
	<p>` + 
          print_cismapper_doc("", 'tracks-bed-short') + `
	</p>
	`;
      break;
  } // switch

  // Return the text or insert it in the element.
  if (id == "") {
    return(html);
  } else {
    document.getElementById(id).insertAdjacentHTML('beforeend', html);
  }
} // print_cismapper_doc

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_cismapper_doc_text(doc_type, extra) {
  var html;
  if (extra == undefined) {extra = ""};

  switch (doc_type) {
    case 'gene-id':
      return(`
	The ID of the gene.
      `);
    case 'gene-name':
      return(`
        The name of the gene. This will be "." if there was no name provided for this 
	gene in your annotation file.
      `);
    case 'tss-id':
      return(`
	The ID of the TSS (transcription start site) of the gene.
        This will be the same as the gene ID if no TSS information was provided in
        your regulatory element (RE) locus file.
      `);
    case 'tss-locus':
      return(`
	The genomic coordinates of the TSS (transcription start site) of the gene.
      `);
    case 'strand':
      return(`
	The chromosomal strand on which the gene is located.
      `);
    case 're-locus':
      return(`
	The genomic coordinates of a (potential) regulatory element (RE) you 
	provided in your (RE) locus file.
      `);
    case 'distance':
      return(`
	The distance between the TSS and the RE_locus, taking Strand into account, 
	so that negative distances mean the RE is upstream of the TSS.
      `);
    case 'histone':
      return(`
	The name of the histone modification used in calculating the
	Pearson correlation between expression of the TSS and
	the level of the histone modification at the RE_Locus.
	The correlation is measured after log-transforming both
	variables: <br>&nbsp;&nbsp;&nbsp;&nbsp; x_new = log(x+1)
      `);
    case 'best-histone':
      return(`
	Equal to 'T' if the link between this RE and TSS based on this histone 
	has the best Score.
      `);
    case 'correlation':
      return(`
	The Pearson correlation between the expression of the 
	TSS and the level of the histone modification at the RE_Locus.
	The correlation is measured after log-transforming both
	variables: <br>&nbsp;&nbsp;&nbsp;&nbsp; x_new = log(x+1)
      `);
    case 'correlation-sign':
      return(`
	The sign of the Pearson correlation of the expression of the 
	TSS and the level of the histone modification at the RE_Locus.
	The correlation is measured after log-transforming both
	variables: <br>&nbsp;&nbsp;&nbsp;&nbsp; x_new = log(x+1)
      `);
    case 'score':
      return(`
	The unadjusted, two-tailed <i>p</i>-value of the Pearson correlation of the
	potential regulatory link between this TSS and RE_Locus.
	The Pearson correlation is between the expression of the 
	TSS and the level of the histone modification at the RE_Locus.
	The correlation is measured after log-transforming both
	variables: <br>&nbsp;&nbsp;&nbsp;&nbsp; x_new = log(x+1)<br>
	The two-tailed <i>p</i>-value is computed assuming the Fisher Z-Transform of
	the correlation coefficient (r) is normally distributed:
	<br>&nbsp;&nbsp;&nbsp;&nbsp; Z = (sqrt(N-3)/2) * (log(1+r) - log(1-r))<br>
	where N is the number of tissues with TSS expression and RE_Locus
	histone levels used in calculating r.
      `);
    case 're-score':
      return(`
	The minimum Score of all potential regulatory links to the TSS,
	adjusted for the number of RE_Locus linked to the TSS.
        The adjustment is <br>&nbsp;&nbsp;&nbsp;&nbsp; new_score = 1 - (1 - score)<sup>n</sup>.
      `);
    case 'best-re':
      return(`
	Equal to 'T' if the link between this RE and TSS has the best RE_Score
	for the TSS.
      `);
    case 'tss-score':
      return(`
	The minimum Score of all potential regulatory links to the TSS,
	adjusted first for the number of RE_Locus linked to the TSS (see RE_Score),
        then for the total number of TSSes (among all Genes) that have links.
        Each adjustment is <br>&nbsp;&nbsp;&nbsp;&nbsp; new_score = 1 - (1 - score)<sup>n</sup>.
      `);
    case 'best-tss':
      return(`
	Equal to 'T' if the link between this RE and TSS has the best TSS_Score
	for the Gene.
      `);
    case 'gene-score':
      return(`
	The minimum Score of all potential regulatory links to the Gene,
	adjusted first for the number of RE_Locus linked to the TSS (see RE_Score), 
	then for the number of TSSes of the Gene, then for the total 
	number of genes that have links.
        Each adjustment is <br>&nbsp;&nbsp;&nbsp;&nbsp; new_score = 1 - (1 - score)<sup>n</sup>.
      `);
    case 'loci-per-tss':
      return(`
	The number of RE_Locus linked to this TSS.
      `);
    case 'num-tsses':
      return(`
	The number of TSSes (of all Genes) that have links.
      `);
    case 'tsses-per-gene':
      return(`
	The number of TSSes of this gene that have links.
      `);
    case 'num-genes':
      return(`
	The number of Genes that have links.
      `);
    case 'gene-table':
      var i = 1;
      return(`
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_ID</td> <td>` + get_cismapper_doc_text('gene-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_Name</td> <td>` + get_cismapper_doc_text('gene-name') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_ID</td> <td>` + get_cismapper_doc_text('tss-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Strand</td> <td>` + get_cismapper_doc_text('strand') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>RE_Locus</td> <td>` + get_cismapper_doc_text('re-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Distance</td> <td>` + get_cismapper_doc_text('distance') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Histone</td> <td>` + get_cismapper_doc_text('histone') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Correlation_Sign</td> <td>` + get_cismapper_doc_text('correlation-sign') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Score</td> <td>` + get_cismapper_doc_text('score') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_Score</td> <td>` + get_cismapper_doc_text('gene-score') + `</td> </tr>
	</table>
      `);
    case 'tss-table':
      var i = 1;
      return(`
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_ID</td> <td>` + get_cismapper_doc_text('tss-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_ID</td> <td>` + get_cismapper_doc_text('gene-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_Name</td> <td>` + get_cismapper_doc_text('gene-name') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Strand</td> <td>` + get_cismapper_doc_text('strand') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>RE_Locus</td> <td>` + get_cismapper_doc_text('re-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Distance</td> <td>` + get_cismapper_doc_text('distance') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Histone</td> <td>` + get_cismapper_doc_text('histone') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Correlation_Sign</td> <td>` + get_cismapper_doc_text('correlation-sign') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Score</td> <td>` + get_cismapper_doc_text('score') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_Score</td> <td>` + get_cismapper_doc_text('tss-score') + `</td> </tr>
	</table>
      `);
    case 'links-table':
      var i = 1;
      return(`
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_ID</td> <td>` + get_cismapper_doc_text('tss-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_Locus</td> <td>` + get_cismapper_doc_text('tss-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Strand</td> <td>` + get_cismapper_doc_text('strand') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>RE_Locus</td> <td>` + get_cismapper_doc_text('re-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Correlation</td> <td>` + get_cismapper_doc_text('correlation') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Score</td> <td>` + get_cismapper_doc_text('score') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>RE_Score</td> <td>` + get_cismapper_doc_text('re-score') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Best_RE</td> <td>` + get_cismapper_doc_text('best-re') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_Score</td> <td>` + get_cismapper_doc_text('tss-score') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Best_TSS</td> <td>` + get_cismapper_doc_text('best-tss') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_Score</td> <td>` + get_cismapper_doc_text('gene-score') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Distance</td> <td>` + get_cismapper_doc_text('distance') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Histone</td> <td>` + get_cismapper_doc_text('histone') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Best_Histone</td> <td>` + get_cismapper_doc_text('best-histone') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_ID</td> <td>` + get_cismapper_doc_text('gene-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_Name</td> <td>` + get_cismapper_doc_text('gene-name') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Loci_per_TSS</td> <td>` + get_cismapper_doc_text('loci-per-tss') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Num_TSSes</td> <td>` + get_cismapper_doc_text('num-tsses') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSSes_per_Gene</td> <td>` + get_cismapper_doc_text('tsses-per-gene') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Num_Genes</td> <td>` + get_cismapper_doc_text('num-genes') + `</td> </tr>
      `);
  } // switch
} // get_cismapper_doc_text

//
// Function to replace the innerHTML of element "id" with an HTML paragraph
// containing the text for 'doc_type', which is known to function get_cismapper_doc_text.
// This function can be used in help pop-ups.
//
function print_cismapper_doc_para(id, doc_type, extra) {
  html = `<p>` + get_cismapper_doc_text(doc_type, extra) + `</p>`;
  document.getElementById(id).insertAdjacentHTML('beforeend', html);
} // print_cismapper_doc_para
