var sort_table = {
  "regulatory_link": [
    {"name": "Gene Score", "fn": sort_gene_score, "priority": 2},
    {"name": "Score", "fn": sort_score, "priority": 1},
    {"name": "Gene ID", "fn": sort_gene_id},
    {"name": "Gene Name", "fn": sort_gene_name},
    {"name": "TSS ID", "fn": sort_tss_id},
    {"name": "RE Locus", "fn": sort_re_locus},
    {"name": "Distance", "fn": sort_distance},
    {"name": "Absolute Distance", "fn": sort_absolute_distance},
    {"name": "Correlation Sign", "fn": sort_correlation_sign},
  ],
};

pre_load_setup();

function pre_load_setup() {
} //pre_load_setup

/*
 * page_loaded
 *
 * Called when the page has loaded for the first time.
 */
function page_loaded() {
  first_load_setup();
  post_load_setup();
} // page_loaded

/*
 * page_loaded
 *
 * Called when a cached page is reshown.
 */
function page_shown(e) {
  if (e.persisted) post_load_setup();
} // page_shown

/*
 * first_load_setup
 *
 * Setup state that is dependent on everything having been loaded already.
 * On browsers which cache state this is only run once.
 */
function first_load_setup() {
    $('show_tss_locus').checked = false;
    $('show_strand').checked = false;
    $('show_histone').checked = false;
    $('show_correlation').checked = false;
    $('show_re_score').checked = false;
    $('show_best_re').checked = false;
    $('show_tss_score').checked = false;
    $('show_best_tss').checked = false;
    $('show_loci_per_tss').checked = false;
    //$('show_num_tsses').checked = false;
    $('show_tsses_per_gene').checked = false;
    //$('show_num_genes').checked = false;
} // first_load_setup

/*
 * post_load_setup
 *
 * Setup state that is dependent on everything having been loaded already.
 */
function post_load_setup() {
  "use strict";
  var tbl, i;

  if (! data['have_tss_info']) {
    toggle_class($("tss_targets_tsv_file"), "hide", 1);
    toggle_class($("tss_elements_tsv_file"), "hide", 1);

    toggle_class($("div_filter_on_best_tss"), "hide", 1);
    toggle_class($("div_filter_on_best_re"), "hide", 1);
    toggle_class($("div_filter_on_tss_id"), "hide", 1);

    toggle_class($("div_show_tss_id"), "hide", 1);
    toggle_class($("div_show_tss_score"), "hide", 1);
    //toggle_class($("div_show_num_tsses"), "hide", 1);
    toggle_class($("div_show_tsses_per_gene"), "hide", 1);
  }

  $("filter_top").disabled = !($("filter_on_top").checked);
  $("filter_best_tss").disabled = !($("filter_on_best_tss").checked);
  $("filter_best_re").disabled = !($("filter_on_best_re").checked);
  $("filter_gene_id").disabled = !($("filter_on_gene_id").checked);
  $("filter_gene_name").disabled = !($("filter_on_gene_name").checked);
  $("filter_tss_id").disabled = !($("filter_on_tss_id").checked);
  $("filter_re_locus").disabled = !($("filter_on_re_locus").checked);
  $("filter_absolute_distance_ge").disabled = !($("filter_on_absolute_distance").checked);
  $("filter_absolute_distance_le").disabled = !($("filter_on_absolute_distance").checked);
  $("filter_distance_ge").disabled = !($("filter_on_distance").checked);
  $("filter_distance_le").disabled = !($("filter_on_distance").checked);
  $("filter_correlation_sign").disabled = !($("filter_on_correlation_sign").checked);
  $("filter_score").disabled = !($("filter_on_score").checked);
  $("filter_gene_score").disabled = !($("filter_on_gene_score").checked);

  tbl = $("regulatory_links");
  toggle_class(tbl, "hide_gene_id", !$("show_gene_id").checked);
  toggle_class(tbl, "hide_gene_name", !$("show_gene_name").checked);
  toggle_class(tbl, "hide_tss_id", !$("show_tss_id").checked || !data["have_tss_info"]);
  toggle_class(tbl, "hide_tss_locus", !$("show_tss_locus").checked);
  toggle_class(tbl, "hide_strand", !$("show_strand").checked);
  toggle_class(tbl, "hide_re_locus", !$("show_re_locus").checked);
  toggle_class(tbl, "hide_distance", !$("show_distance").checked);
  toggle_class(tbl, "hide_histone", !$("show_histone").checked);
  toggle_class(tbl, "hide_correlation", !$("show_correlation").checked);
  toggle_class(tbl, "hide_correlation_sign", !$("show_correlation_sign").checked);
  toggle_class(tbl, "hide_score", !$("show_score").checked);
  toggle_class(tbl, "hide_re_score", !$("show_re_score").checked);
  toggle_class(tbl, "hide_best_re", !$("show_best_re").checked);
  toggle_class(tbl, "hide_tss_score", !$("show_tss_score").checked || !data["have_tss_info"]);
  toggle_class(tbl, "hide_best_tss", !$("show_best_tss").checked);
  toggle_class(tbl, "hide_gene_score", !$("show_gene_score").checked);
  toggle_class(tbl, "hide_loci_per_tss", !$("show_loci_per_tss").checked);
  //toggle_class(tbl, "hide_num_tsses", !$("show_num_tsses").checked || !data["have_tss_info"]);
  toggle_class(tbl, "hide_tsses_per_gene", !$("show_tsses_per_gene").checked || !data["have_tss_info"]);
  //toggle_class(tbl, "hide_num_genes", !$("show_num_genes").checked);

  make_regulatory_links_table();
} // post_load_setup

/*
 * toggle_filter
 *
 * Called when the user clicks a checkbox
 * to enable/disable a filter option.
 */
function toggle_filter(chkbox, filter_id) {
  var filter = $(filter_id);
  filter.disabled = !(chkbox.checked);
  if (!filter.disabled) {
    filter.focus();
    if (filter.select) filter.select();
  }
} // toggle_filter

/*
 * enable_filter
 *
 * Called when the user clicks a filter label.
 * Enables the filter.
 */
function enable_filter(chkbox_id, filter_id) {
  var chkbox = $(chkbox_id);
  if (!chkbox.checked) {
    var filter = $(filter_id);
    $(chkbox_id).checked = true;
    filter.disabled = false;
    filter.focus();
    if (filter.select) filter.select();
  }
} // enable_filter

/*
 * update_filter
 *
 * If the key event is an enter key press then
 * update the filter on the regulatory links table
 */
function update_filter(e) {
  if (!e) var e = window.event;
  var code = (e.keyCode ? e.keyCode : e.which);
  if (code == 13) {
    e.preventDefault();
    make_regulatory_links_table();
  }
} // update_filter

function num_keys(e) {
  if (!e) var e = window.event;
  var code = (e.keyCode ? e.keyCode : e.which);
  var keychar = String.fromCharCode(code);
  var numre = /\d/;
  // only allow 0-9 and various control characters (Enter, backspace, delete)
  if (code != 8 && code != 9 && code != 13 && code != 46 && !numre.test(keychar)) {
    e.preventDefault();
  }
}

/*
 * toggle_column
 *
 * Adds or removes a class from the table displaying the
 * predicted regulatory links. This is primary used to set the visibility
 * of columns by using css rules. If the parameter 'show' is not passed
 * then the existence of the class will be toggled, otherwise it will be
 * included if show is false.
 */
function toggle_column(cls) {
  toggle_class($("regulatory_links"), cls);
} // toggle_column

function populate_sort_list (sellist, items) {
  var i, j, item, opt, priority, selected;
  priority = 0;
  selected = 0;
  for (i = 0, j = 0; i < items.length; i++) {
    item = items[i];
    if (typeof item["show"] === 'undefined' || item["show"]) {
      opt = document.createElement("option");
      opt.innerHTML = item["name"];
      opt.value = i;
      sellist.add(opt, null);
      if (typeof item["priority"] !== 'undefined' && item["priority"] > priority) {
        selected = j;
        priority = item["priority"];
      }
      j++;
    }
  }
  sellist.selectedIndex = selected;
} // populate_sort_list

function populate_sort_lists() {
  "use strict";
  var i, motif_sort, regulatory_link_sort_list;
  regulatory_link_sort = $("regulatory_link_sort");
  //try {regulatory_link_sort.removeEventListener("click", sync_sort_selection, false);} catch (e) {}
  for (i = regulatory_link_sort.options.length-1; i >= 0; i--) regulatory_link_sort.remove(i);
  regulatory_link_sort_list = sort_table["regulatory_link"];
  populate_sort_list(regulatory_link_sort, regulatory_link_sort_list);
  //sync_sort_selection();
  //regulatory_link_sort.addEventListener("click", sync_sort_selection, false);
  //regulatory_link_sort.addEventListener("change", sync_sort_selection, false);
} // populate sort lists

// Sort by increasing Gene ID, then by TSS ID, then RE.
function sort_gene_id(link1, link2) {
  var diff;
  diff = link1['gene_id'].localeCompare(link2['gene_id']);
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_gene_id

// Sort by increasing Gene Name, then by TSS ID, then RE.
function sort_gene_name(link1, link2) {
  var diff;
  diff = link1['gene_name'].localeCompare(link2['gene_name']);
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_gene_name

// Sort by increasing TSS ID, then RE.
function sort_tss_id(link1, link2) {
  var diff;
  diff = link1['tss_id'].localeCompare(link2['tss_id']);
  if (diff == 0) { 
    diff = link1['re_locus'].localeCompare(link2['re_locus']);
  }
  return diff;
} // sort_tss_id

// Sort by increasing RE, then by TSS ID.
function sort_re_locus(link1, link2) {
  var diff;
  diff = link1['re_locus'].localeCompare(link2['re_locus']);
  if (diff == 0) { 
    diff = link1['tss_id'].localeCompare(link2['tss_id']);
  }
  return diff;
} // 

// Sort by increasing distance, then TSS, then RE.
function sort_distance(link1, link2) {
  var diff;
  diff = link1['distance'] - link2['distance'];
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_distance

// Sort by increasing absolute distance, then TSS, then RE.
function sort_absolute_distance(link1, link2) {
  var diff;
  diff = Math.abs(link1['distance']) - Math.abs(link2['distance']);
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_absolute_distance

// Sort by correlation sign ("+" comes first), then TSS, then RE.
function sort_correlation_sign(link1, link2) {
  var diff, sign1, sign2;
  sign1 = link1['correlation'] > 0 ? '+' : '-';
  sign2 = link2['correlation'] > 0 ? '+' : '-';
  diff = sign1.localeCompare(sign2);
  if (diff != 0) return -diff;
  return sort_tss_id(link1, link2);
} // sort_correlation_sign

// Sort by increasing gene_score, then score, then TSS, then RE.
function sort_gene_score(link1, link2) {
  var diff;
  diff = link1['gene_score'] - link2['gene_score'];
  if (diff != 0) return diff;
  diff = link1['score'] - link2['score'];
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_gene_score

// Sort by increasing score, then TSS, then RE.
function sort_score(link1, link2) {
  var diff;
  diff = link1['score'] - link2['score'];
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_score

/*
 * regulatory_link_sort_cmp
 *
 * Gets the sorting comparator by index.
 *
 */
function regulatory_link_sort_cmp(index) {
  "use strict";
  if (index > 0 && index < sort_table["regulatory_link"].length) {
    return sort_table["regulatory_link"][index]["fn"];
  }
  return sort_gene_score;	// default sorting function
}

function get_filter() {
  var filter, pat, value, value1, value2, count;

  filter = {};

  // get the regulatory link count limit
  filter["on_count"] = $("filter_on_top").checked;
  count = parseFloat($("filter_top").value);
  if (isNaN(count) || count < 1) {
    filter["on_count"] = false;
    $("filter_top").className = "error";
  } else {
    filter["count"] = count;
    $("filter_top").className = "";
  }

  // get the best TSS filter
  filter["on_best_tss"] = $("filter_on_best_tss").checked;

  // get the best RE filter
  filter["on_best_re"] = $("filter_on_best_re").checked;

  // get the gene id filter
  filter["on_gene_id"] = $("filter_on_gene_id").checked;
  pat = $("filter_gene_id").value;
  try {
    filter["gene_id"] = new RegExp(pat, "i");
    $("filter_gene_id").className = "";
  } catch (err) {
    $("filter_gene_id").className = "error";
    filter["on_gene_id"] = false;
  }

  // get the gene name filter
  filter["on_gene_name"] = $("filter_on_gene_name").checked;
  pat = $("filter_gene_name").value;
  try {
    filter["gene_name"] = new RegExp(pat, "i");
    $("filter_gene_name").className = "";
  } catch (err) {
    filter["on_gene_name"] = false;
    $("filter_gene_name").className = "error";
  }

  // get the TSS id filter
  filter["on_tss_id"] = $("filter_on_tss_id").checked;
  pat = $("filter_tss_id").value;
  try {
    filter["tss_id"] = new RegExp(pat, "i");
    $("filter_tss_id").className = "";
  } catch (err) {
    $("filter_tss_id").className = "error";
    filter["on_tss_id"] = false;
  }

  // get the RE locus filter
  filter["on_re_locus"] = $("filter_on_re_locus").checked;
  pat = $("filter_re_locus").value;
  try {
    filter["re_locus"] = new RegExp(pat, "i");
    $("filter_re_locus").className = "";
  } catch (err) {
    $("filter_re_locus").className = "error";
    filter["on_re_locus"] = false;
  }

  // get the absolute distance filter
  filter["on_absolute_distance"] = $("filter_on_absolute_distance").checked;
  value1 = parseFloat($("filter_absolute_distance_ge").value);
  value2 = parseFloat($("filter_absolute_distance_le").value);
  if (isNaN(value1)) {
    $("filter_absolute_distance_ge").className = "error";
  } else if (isNaN(value2)) {
    $("filter_absolute_distance_le").className = "error";
  } else {
    if (value1 < 0 || value1 > value2) {
      $("filter_absolute_distance_ge").className = "error";
    }
    if (value2 < 0) {
      $("filter_absolute_distance_le").className = "error";
    }
  }
  if (isNaN(value1) || isNaN(value2) || value1 < 0 || value2 < 0 || value1 > value2) {
    filter["on_absolute_distance"] = false;
    filter["on_absolute_distance_ge"] = false;
    filter["on_absolute_distance_le"] = false;
  } else {
    filter["on_absolute_distance_ge"] = true;
    $("filter_absolute_distance_ge").className = "";
    filter["on_absolute_distance_le"] = true;
    $("filter_absolute_distance_le").className = "";
    filter["absolute_distance_ge"] = value1;
    filter["absolute_distance_le"] = value2;
  }

  // get the distance filter
  filter["on_distance"] = $("filter_on_distance").checked;
  value1 = parseFloat($("filter_distance_ge").value);
  value2 = parseFloat($("filter_distance_le").value);
  if (isNaN(value1)) {
    $("filter_distance_ge").className = "error";
  } else if (isNaN(value2)) {
    $("filter_distance_le").className = "error";
  } else {
    if (value1 > value2) {
      $("filter_distance_ge").className = "error";
    }
  }
  if (isNaN(value1) || isNaN(value2) || value1 > value2) {
    filter["on_distance"] = false;
    filter["on_distance_ge"] = false;
    filter["on_distance_le"] = false;
  } else {
    filter["on_distance_ge"] = true;
    $("filter_distance_ge").className = "";
    filter["on_distance_le"] = true;
    $("filter_distance_le").className = "";
    filter["distance_ge"] = value1;
    filter["distance_le"] = value2;
  }

  // get the correlation sign filter
  filter["on_correlation_sign"] = $("filter_on_correlation_sign").checked;
  value = $("filter_correlation_sign").value;
  if (value == "") {
    $("filter_correlation_sign").className = "error";
    filter["on_correlation_sign"] = false;
  } else {
    filter["correlation_sign"] = value[0];
    $("filter_correlation_sign").className = "";
  }

  // get the score filter
  filter["on_score"] = $("filter_on_score").checked;
  if ((value = parseFloat($("filter_score").value)) != null) {
    filter["score"] = value;
    $("filter_score").className = "";
  } else {
    filter["score"] = false;
    $("filter_score").className = "error";
  }

  // get the gene score filter
  filter["on_gene_score"] = $("filter_on_gene_score").checked;
  if ((value = parseFloat($("filter_gene_score").value)) != null) {
    filter["gene_score"] = value;
    $("filter_gene_score").className = "";
  } else {
    filter["gene_score"] = false;
    $("filter_gene_score").className = "error";
  }

  return filter;
} //get_filter

function filter_regulatory_link(filter, regulatory_link) {
  if (filter["on_best_tss"] && (regulatory_link["best_tss"]) != "T") return true;
  if (filter["on_best_re"] && (regulatory_link["best_re"]) != "T") return true;
  if (filter["on_gene_id"] && !filter["gene_id"].test(regulatory_link["gene_id"])) return true;
  if (filter["on_gene_name"] && !filter["gene_name"].test(regulatory_link["gene_name"])) return true;
  if (filter["on_tss_id"] && !filter["tss_id"].test(regulatory_link["tss_id"])) return true;
  if (filter["on_re_locus"] && !filter["re_locus"].test(regulatory_link["re_locus"])) return true;
  if (filter["on_absolute_distance"] && (
    (Math.abs(regulatory_link["distance"]) < filter["absolute_distance_ge"]) ||
    (Math.abs(regulatory_link["distance"]) > filter["absolute_distance_le"]))
  ) return true;
  if (filter["on_distance"] && (
    (regulatory_link["distance"] < filter["distance_ge"]) ||
    (regulatory_link["distance"] > filter["distance_le"]))
  ) return true;
  if (filter["on_correlation_sign"] && (regulatory_link["correlation_sign"] != filter["correlation_sign"])) return true;
  if (filter["on_score"] && regulatory_link["score"] > filter["score"]) return true;
  if (filter["on_gene_score"] && regulatory_link["gene_score"] > filter["gene_score"]) return true;
  return false;
}

function filter_regulatory_links_table() {
  "use strict";
  var tbl, tbody, tr, regulatory_link, i, filter;
  tbl = $("regulatory_links");
  filter = get_filter();
  for (i = 0; i < tbl.tBodies.length; i++) {
    tbody = tbl.tBodies[i];
    regulatory_link = tbody["data_regulatory_link"];
    toggle_class(tbody, "filtered", filter_motif(filter, regulatory_link));
  }
}

function make_other_settings() {
  $("opt_rna_source").textContent = data.options.rna_source;
  $("opt_tissues").textContent = data.options.tissues;
  $("opt_histone_root").textContent = data.options.histone_root;
  $("opt_histone_names").textContent = data.options.histone_names;
  $("opt_max_link_distances").textContent = data.options.max_link_distances;
  $("opt_expression_root").textContent = data.options.expression_root;
  $("opt_expression_file_type").textContent = data.options.expression_file_type;
  $("opt_transcript_types").textContent = data.options.transcript_types;
  $("opt_min_feature_count").textContent = data.options.min_feature_count;
  $("opt_min_max_expression").textContent = data.options.min_max_expression;
  $("opt_max_html_score").textContent = data.options.max_html_score;
} // make_other_settings

function make_regulatory_links_row(tbody, link) {
  var row;
  row = tbody.insertRow(tbody.rows.length);
  add_text_cell(row, link['gene_id'], 'col_gene_id');
  add_text_cell(row, link['gene_name'], 'col_gene_name');
  add_text_cell(row, link['tss_id'], 'col_tss_id');
  add_text_cell(row, link['tss_locus'], 'col_tss_locus');
  add_text_cell(row, link['strand'], 'col_strand');
  add_text_cell(row, link['re_locus'], 'col_re_locus');
  add_text_cell(row, link['distance'], 'col_distance');
  add_text_cell(row, link['histone'], 'col_histone');
  add_text_cell(row, link['correlation'].toFixed(5), 'col_correlation');
  add_text_cell(row, link['correlation'] > 0 ? '+' : '-', 'col_correlation_sign');
  add_text_cell(row, link['score'].toExponential(3), 'col_score');
  add_text_cell(row, link['re_score'].toExponential(3), 'col_re_score');
  add_text_cell(row, link['best_re'], 'col_best_re');
  add_text_cell(row, link['tss_score'].toExponential(3), 'col_tss_score');
  add_text_cell(row, link['best_tss'], 'col_best_tss');
  add_text_cell(row, link['gene_score'].toExponential(3), 'col_gene_score');
  add_text_cell(row, link['loci_per_tss'], 'col_loci_per_tss');
  //add_text_cell(row, link['num_tsses'], 'col_num_tsses');
  add_text_cell(row, link['tsses_per_gene'], 'col_tsses_per_gene');
  //add_text_cell(row, link['num_genes'], 'col_num_genes');
} // make_regulatory_links_row

function make_regulatory_links_table() {
  var regulatory_link_sort;
  var filter, regulatory_links, filtered, regulatory_link;
  var tbl, tbody, row, cell, i, skipped;

  // get the filter and sort comparator
  filter = get_filter();
  regulatory_link_sort = regulatory_link_sort_cmp(parseInt($('regulatory_link_sort').value));
 
  // filter the regulatory links
  regulatory_links = data['regulatory_links']; 
  filtered = [];
  for (i = 0; i < regulatory_links.length; i++) {
    if (filter_regulatory_link(filter, regulatory_links[i])) continue;
    if (regulatory_links[i]['filtered']) continue;
    filtered.push(regulatory_links[i]);
  }
  filtered.sort(regulatory_link_sort);
  // limit
  if (filter['on_count']) {
    if (filtered.length > filter['count']) filtered.length = filter['count'];
  }

  // clear the table
  tbl = $('regulatory_links');
  for (i = tbl.tBodies.length - 1; i >= 0; i--) {
    tbody = tbl.tBodies[i];
    tbody.parentNode.removeChild(tbody);
  }

  // add the new rows to the table
  for (i = 0; i < filtered.length; i++) {
    regulatory_link = filtered[i];
    tbody = document.createElement('tbody');
    tbody.className = "regulatory_link_group";
    tbody['data_regulatory_link'] = regulatory_link;
    make_regulatory_links_row(tbody, regulatory_link);
    tbl.appendChild(tbody);
  }

  // note the count of filtered regulatory links
  if (filtered.length != regulatory_links.length) {
    skipped =  regulatory_links.length - filtered.length;
    tbody = document.createElement('tbody');
    row = tbody.insertRow(tbody.rows.length);

    if (skipped === 1) {
      desc = "1 regulatory link hidden due to filters";
    } else {
      desc = skipped + " regulatory links hidden due to filters";
    }
    cell = row.insertCell(row.cells.length);
    cell.colSpan = 19;
    cell.style.textAlign = "center";
    cell.style.fontWeight = "bold";
    cell.appendChild(document.createTextNode(desc));
    tbl.appendChild(tbody);
  }
} // make_regulatory_links_table
