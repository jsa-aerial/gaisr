
// Set both the field and cookie.  Generally this should only be
// called once per new user site useage.  This is because our init on
// site page load checks for the user cookie and sets the field
//
function setUser (u) {
    $('user').setValue(u);
    setCookie("user", u, 365);
}


// Setup a basic GreyBox window.  TTITLE is the title of the GB
// window, W/H give the basic size of the window, DIV is the main
// container of window, divID is the element ID of scaffholding div of
// the window content, and divContent is the content element put in
// divID element.
//
function setupGBox (title, w, h, div, divID, divContent) {
    GB_showCenter(title, '', h, w);
    var gbw =  $('GB_window');
    gbw.children[1].update(div);
    var gbdiv = $(divID);
    gbdiv.update(divContent);
}




/*
  -- Page tabbing
 */

var pageInfo = {
    1: {elts: ["seq-query", "actionForm"],
        results: "", cnt: 0, tbox: false},
    2: {elts: ["result-view", "evspan", "view-div"],
        results: "", cnt: 0, tbox: false},
    3: {elts: [], results: "", cnt: 0},
    4: {elts: [], results: "", cnt: 0}
}

function setDisplay (o, n) {
    o.each(function(e){$(e).hide()});
    n.each(function(e){$(e).show()});
}

function recallResults(tp) {
    var n = pageInfo[tp];
    var results = $('results');
    results.update("");
    $('toggleBox').setValue(n.tbox);
    if (n.results == "") {
        $('count').hide();
    } else {
        n.results.each(function (i) {results.appendChild(i);});
        var count = $('count');
        count.update(n.cnt);
        count.show();
    }
}

var activeTab = undefined;

function setActiveTab (elt) {
    var tmp = activeTab;
    tmp.removeClassName('tab-active');
    tmp.addClassName('tab-inactive');
    elt.removeClassName('tab-inactive');
    elt.addClassName('tab-active');
    activeTab = elt;

    var tabPage = elt.getAttribute('value');
    var oldPage = tmp.getAttribute('value');
    pageInfo[oldPage].results = $('results').childElements();
    pageInfo[oldPage].cnt = $('count').textContent;
    pageInfo[oldPage].tbox = $('toggleBox').getValue();
    setDisplay(pageInfo[oldPage].elts, pageInfo[tabPage].elts);
    recallResults(tabPage);
}




/*
  --- Primary Query Request Processing ---
*/


function createMarkerElement(item) {
    var info = {name: item.taxname, taxid: item.taxon_id}
    var elt = $(document.createElement('li'));
    elt.update(Templates['marker'].evaluate(info));
    return elt;
}


function createResultElement (item, template) {
    var elt = $(document.createElement('li'));
    elt.update(Templates[template].evaluate(item));
    return elt;
}


function itemCompare (l, r) {
    // This is a completely idiotic and bogus aspect of JS that it
    // requires this sort of nonsense
    var ll = l.taxname.toLowerCase();
    var lr = r.taxname.toLowerCase();
    if (ll < lr) {
        return -1;
    } else if (ll == lr) {
        return 0;
    } else {
        return 1;
    }
}


function populateResults(event) {
    //console.log("Start Populate Results");
    var template = event.memo.template;
    var items = event.memo.items.sort(itemCompare);
    var results = $('results');
    var marker = "";
    var markerID = "";
    items.each(function (item, index) {
        if (marker != item.taxname) {
            marker = item.taxname;
            markerID = item.taxon_id;
            var mElt = createMarkerElement(item);
            results.appendChild(mElt);
        }
        var elt = createResultElement(item, template);
        results.appendChild(elt);
    });
    Accordion.init($$('#results'));
    //console.log("***DONE");
}


function createIMPair (item) {
    var pair = $(document.createElement('li'));
    pair.update(Templates['chat'].evaluate(item));
    return pair;
}

function processResults (request) {
    //console.log("Start Process Results");
    var items = request.responseJSON;
    var results = $('results');
    var chat = $('chat');
    var size = {count: 0};
    size.count = items[0];
    $('toggleBox').setValue(false);
    $('count').update(Templates['count'].evaluate(size));
    chat.insertBefore(createIMPair(items[1]),chat.down());
    $('query').clear();
    results.update("");
    document.fire("results:updated",
                  {items: items.slice(2), template: 'qentry'});
    //console.log("***DONE");
}




function makeQuery (form, e) {
    var u = $('user').getValue();

    if (u == "") {
        u = prompt('Query is missing USER', "enter user name here");
        if (u == "") {
            alert('Query requires USER');
            return undefined;
        }
        setUser(u);
    }

    var user = form.children[0];
    user.setValue(u);

    // Update input history
    addHistoryItem($('query').getValue());

    // Actually, we don't need sleep here, but for consistent form of
    // how we are firing off Ajax requests, set things up like for all
    // the others and fire it off with setTimeout...
    //
    var launchAjax = function () {
        new Ajax.Request(form.action, {
            method: 'get',
            parameters: form.serialize(),
            onSuccess: function (request) {
                processLog(request);
                processResults(request);
            }
        });
    }

    setTimeout(launchAjax, 10); // Effectively "instantaneous"...
}




/*
  --- Feature fetch and display ---
*/


function createFeatureElement(item) {
    var elt = $(document.createElement('li'));
    elt.update(Templates['feature'].evaluate(item));
    //console.log("createFeatureElement, item = " + item + ", elt = " + elt);
    return elt;
}


var itemDbg = "NA";

function addLocsNVS (item) {
    //console.log("Start addLocsNVS");
    itemDbg = item;
    var sfid = item.sfid;
    var sftype = item.sftype;
    var locs = item.locs;
    var nvs = item.nvs;
    var loc_ul = $(sfid+"loc");
    var nvs_ul = $(sfid+"nvs");
    locs.each(function (loc) {
        var elt = $(document.createElement('li'));
        elt.update(Templates['loc'].evaluate(loc));
        loc_ul.appendChild(elt);
    });
    nvs.each(function (nv) {
        var x = nv;
        var elt = $(document.createElement('li'));
        if (x.name == "translation")
            x.value = x.value.substring(0, [x.value.length, 200].min());
        elt.update(Templates['nv'].evaluate(x));
        nvs_ul.appendChild(elt);
    });
    //console.log("***DONE");
}


function processFeatures (request) {
    //console.log("Start processFeatures");
    reqJson = request;
    var size = {count: 0};
    var items = request.responseJSON;
    var gbid = items[0];
    var fresults = $('fresults');

    items = items.slice(1);
    size.count = items.length;
    $('fcount').update(Templates['count'].evaluate(size));

    items.each(function (item, index) {
        var elt = createFeatureElement(item);
        fresults.appendChild(elt);
        addLocsNVS(item);
    });
    Accordion.init($$('#fresults'));
    //console.log("***DONE");
}


function createFeatureGBWindow (title, fname, fid, gbid) {
    listGBox(title, 800, 700);
    var launchAjax = function () {
        new Ajax.Request("/mlab/features", {
            method: 'get',
            parameters: {id: fid, name: fname, gbid: gbid},
            onSuccess: function (request) {
                //processLog(request)
                processFeatures(request);
            }
        });
    }
    // This looks crazy, but it effectively acts as a sleep which
    // gives the greybox time to fully paint before Ajax response
    // takes over
    setTimeout(launchAjax, 800);
}




var scriblInfo = {
    width: 0,
    tipStyle: "light",
    featHeight: 23,
    featTextSize: 12,
    featColors: {hit: "red",
                 gene: "rgb(163, 127, 163)", // "rgb(123, 107, 143)"
                 CDS: "rgb(233, 151, 63)",
                 misc_feature: "rgb(203, 95, 93)"},

    featNames: {hit: "loc", gene: "gene", CDS: "gene", misc_feature: "loc"},

    featNotes: {hit: ["locus"],
                gene: ["gene", "locus", "note"],
                CDS:  ["gene", "locus", "product", "protein_id", "note"],
                misc_feature: ["locus", "note"]},

    featSeq: [],
    chartInfos: {}
}


function featureURL (gbid, start, end) {
    var urlBase = "http://www.ncbi.nlm.nih.gov/nuccore/";
    var gbid = gbid + "?";
    var startArg = "from=" + start + "&";
    var endArg = "to=" + end + "&";
    return urlBase + gbid + startArg + endArg + "report=gbwithparts";
}


function scriblLocusInfo (name, item, chart) {
    //console.log("Start ScribleLocusInfo");
    var chartInfo = scriblInfo.chartInfos[name];

    var sfid = item.sfid;
    var sftype = item.sftype;
    var locs = item.locs;

    var gbid = chartInfo.gbid;
    var scriblCkBoxMap = chartInfo.scriblCkBoxMap;
    var colors = scriblInfo.featColors;
    var names  = scriblInfo.featNames;
    var noteNames  = scriblInfo.featNotes;

    var nvs = item.nvs;
    var nvm = {}
    nvs.each(function(nv) {nvm[nv.name]=nv.value;});

    var getVal = function (k){
        var v = nvm[k] || "";
        var n = (k == "gene") ? "" : k + ": ";
        return (v != "") ? n + v : "";
    }

    var locSeq = [];
    locs.each(function (loc) {
        var box = scriblCkBoxMap[sftype];
        if (sftype == "hit" || (box != undefined && box.getValue() == "on")) {
            var start = loc.start;
            var len = loc.len;
            var strand = (loc.strand == 1) ? '+' : '-';
            //console.log("loc.strand: " + loc.strand + ", strand: " + strand);
            locSeq.push(chart.addGene(start, len, strand));

            var feature = locSeq[locSeq.length-1];
            feature.color = colors[sftype];

            feature.name = names[sftype];
            if (feature.name != "loc") {
                feature.name = nvm[feature.name] || "NA";
            } else {
                feature.name = start+".."+loc.end
            }

            var notes = noteNames[sftype].map(
                function(x){
                    return (x == "locus") ? start+".."+loc.end : getVal(x)});
            //console.log("Notes: " + notes.join(", "));
            notes = notes.without("").join(", ");
            feature.onMouseover = notes;
            feature.onClick = featureURL(gbid, start, loc.end);
        }
    });
    scriblInfo.featSeq = locSeq;
    //console.log("***DONE scriblLocusInfo");
}


function scriblFeatureTypes (name, info) {
    //console.log("Start ScriblFeatureTypes");
    var w = info.width;
    var chartInfo = info.chartInfos[name];
    var items = chartInfo.featureSet;
    var canvas = chartInfo.canvas;
    var chart = new Scribl(canvas, w - 20);

    if (chartInfo.chart) chartInfo.chart.tracks = []; // Clear old chart!
    chartInfo.chart = chart;
    chart.tooltips.style = scriblInfo.tipStyle;
    chart.trackSizes = scriblInfo.featHeight;
    chart.gene.text.size = scriblInfo.featTextSize;
    items.each(function (item) {
            scriblLocusInfo(name, item, chart);
    });
    chart.redraw();
    //console.log("***DONE scriblFeatureTypes");
}


function reScribl (elt) {
    var name = elt.getAttribute('name');
    scriblFeatureTypes(name, scriblInfo);
}


function lociCount (items) {
    var count = 0;
    items.each(function (item) {
            count = count + item.locs.length;
        });
    return count;
}

function scriblFeatures (name, request, canvas) {
    //console.log("Start Scrible Features");
    scriblInfo.chartInfos[name] = {}
    var chartInfo = scriblInfo.chartInfos[name];
    chartInfo.canvas = canvas;

    var w = canvas.getWidth();
    var h = canvas.getHeight();
    var size = {feats: 0, loci: 0};
    var items = request.responseJSON;
    scriblInfo.width = w;

    size.feats = items.length;
    size.loci = lociCount(items);
    $(name+'fcount').update(Templates['scribcount'].evaluate(size));
    $(name+'fcount').show();

    chartInfo.gbid = request.gbid || "";
    chartInfo.featureSet = items;
    chartInfo.scriblCkBoxMap = {gene: $(name+'gbx'),
                                CDS: $(name+'cbx'),
                                misc_feature: $(name+'mbx')};
    $(name+'gbx').setValue(true);
    scriblFeatureTypes(name, scriblInfo);
    //console.log("***DONE");
}


function createScriblGBWindow (title, fname, fid, gbid) {
    var gbW = 800;
    var gbH = 600;
    var scriblW = gbW - 20;
    var scriblH = gbH - 20;
    scriblGBox(title, gbW, gbH);

    var launchAjax = function () {
        scriblInfo.charInfos = {};
        var ul = $('chartUL');
        var elt = $(document.createElement('li'));
        var scriblElt = Templates['scribl'].evaluate(
            {name: fname, hitloc: "na", scriblW: scriblW, scriblH: scriblH});
        elt.update(scriblElt);
        ul.appendChild(elt);
        var canvas = $(fname); // It's now a DOM element!
        new Ajax.Request("/mlab/features", {
            method: 'get',
            parameters: {id: fid, name: fname, gbid: gbid},
            onSuccess: function (request) {
                processLog(request);
                // First element is the GID for use in NCBI queries
                request.gbid = request.responseJSON[0];
                request.responseJSON = request.responseJSON.slice(1);
                scriblFeatures(fname, request, canvas);
            }
        });
    }
    // This looks crazy, but it effectively acts as a sleep which
    // gives the greybox time to fully paint before Ajax response
    // takes over
    setTimeout(launchAjax, 800);
}




// View/process Hit File Uploads --------------------------------------------

function scriblToggleChecks (mbx) {
    console.log("Mark Genomes checked -- " + mbx.getValue());
    var boxes = $$('.scrib-item-box');
    if (mbx.getValue() == null) {
        boxes.each(function(b){b.setValue(false);});
    } else {
        boxes.each(function(b){b.setValue(true);});
    }
}


function startUpload() {
    $('count').hide();
    $('evspan').hide();
    $('spinner').show();
}


var hitFeatures = {};

function uploadDone() { //Function will be called when hidden-iframe loaded
    $('toggleBox').setValue(false);
    $('spinner').hide();

    var respIFrame = frames['hidden-iframe'];
    var resJSON = respIFrame.document.getElementsByTagName("body")[0].innerHTML;
    var items = eval("(" + resJSON + ")");

    if (items.error) {
        alert(items.error);
        return undefined;
    } else if (items.info) {
        alert(items.info);
        return undefined;
    }

    var size = {count: items.length};
    $('count').update(Templates['count'].evaluate(size));
    $('count').show();
    $('evspan').show();
    $('results').update("");
    hitFeatures = {};
    items.each(function(entry){
            hitFeatures[entry.name] = {hitloc: entry.hitloc,
                                       hitstrand: entry.hitstrand,
                                       hitfeats: entry.hit_features,
                                       evalue: parseFloat(entry.evalue),
                                       gbid: entry.gbid,
                                       taxname: entry.taxname,
                                       ancestors: entry.ancestors};
        });
    populateResults({memo: {items: items, template: 'ventry'}});
}


function evalueCompare (l, r) {
    var el = hitFeatures[l]["evalue"];
    var er = hitFeatures[r]["evalue"];
    if (el < er) {return -1;} else if (el == er) {return 0;} else {return 1;}
}


function scriblRequest (e) {
    console.log("scriblRequest: " + e);
    var gbW = 940;
    var gbH = 816;
    var scriblW = gbW - 40;
    var scriblH = [gbH / 3, 300].min();
    var entries = gatherChecked('.item-box').sort(evalueCompare);

    scriblGBox(entries.slice(0,7).join(","), gbW, gbH);
    scriblInfo.chartInfos = {};

    var count = 0
    var nlfmap = entries.map(function(e){
            var m = {}
            m.name = e;
            m.loc = hitFeatures[e]["hitloc"]
            m.features = hitFeatures[e]["hitfeats"];
            m.evalue = hitFeatures[e]["evalue"];
            m.gbid = hitFeatures[e]["gbid"];
            count = count + m.features.length;
            return m;
        });
    console.log("COUNT: " + count);

    var scriblEntries = function () {
        $('scriblCtrl').show();
        var ul = $('chartUL');
        nlfmap.each(function(nlf){
            // First set event watchers
            $('genFastaForm').observe('submit', function(e){
                    stopEvent(e)});
            $('nameSaveForm').observe('submit', function(e){
                    stopEvent(e);});
            // Fill in and paint canvas
            var fname = nlf.name;
            var request = {responseJSON: nlf.features, gbid: nlf.gbid};
            var elt = $(document.createElement('li'));
            var scriblElt = Templates['scribl'].evaluate(
                {name: fname, hitloc: nlf.loc, evalue: nlf.evalue,
                 scriblW: scriblW, scriblH: scriblH});
            elt.update(scriblElt);
            ul.appendChild(elt);
            var canvas = $(fname); // It's now a DOM element!
            scriblFeatures(fname, request, canvas);
            //$(fname+'fcount').show();
            });
    }

    // This is INSANE, but if we don't hack this into a timed request,
    // it doesn't work.  I think it has to do with some sort of timing
    // of rendering the GBox and its DOM.  JavaScript - weep for all
    // who suffer it!!
    setTimeout(scriblEntries, 800);
}


function genomeNamePage () {
    console.log("genomeNamePage ...");
    var gndiv = $('genomeNamesDiv');
    if (gndiv.getAttribute('style') == "") {
        gndiv.hide();
    } else {
        var entries = gatherChecked('.scrib-item-box');
        entries = entries.map(function(e) {
                var name = e.slice(0, e.search("Scrib"));
                return name + ": " + hitFeatures[name]["taxname"];
            });
        $('genomeNames').setValue(entries.join('\n'));
        $('genomeNamesDiv').show();
        //window.open(entries.join('\n'),"genomes");
    }
}


function scriblGenFasta() {
    var form = $('scriblForm');
    var entries = gatherChecked('.scrib-item-box');
    var s = entries.map(function(e) {
            var name = e.slice(0, e.search("Scrib"));
            return name + " " + hitFeatures[name]["hitloc"] +
                   "/" + hitFeatures[name]["hitstrand"] +
                   ": ev(" + hitFeatures[name]["evalue"] + "), " +
                   hitFeatures[name]["ancestors"];
        });
    var u = $('user').getValue();
    var n = $('sfname').getValue();

    if (u == "" || n == "") {
        alert('Request requires both USER and Filename');
        if (u == "") {
            u = prompt("Enter user", "enter user name here");
            if (u == "") {
                alert(a + 'requires user');
                return undefined;
            }
            setUser(u);
        }
        if (n == "") return undefined;
    }

    var user = form.children[0];
    var act = form.children[1];
    var txt = form.children[2];
    var fname = form.children[3];
    user.setValue(u);
    act.setValue('genfasta');
    txt.setValue(s.join("$$")); // Can't use new line as Chrome,FF4 -> space!
    fname.setValue(n);

    var launchAjax = function () {
        new Ajax.Request(form.action, {
            method: 'post',
            parameters: form.serialize(),
            onSuccess: function (request) {
                processLog(request);
                alert("Generated:\n" + request.responseJSON);
            }
        });
    }
    setTimeout(launchAjax, 500);
}


function dbSaveNames () {
    alert('Saving selected item names to DB Not Yet Implemented');
}



var imgPNGs = undefined;

function genGenomePNG () {
    console.log("genGenomePNG ...");
    var entries = gatherChecked('.scrib-item-box');
    entries = entries.map(function(e) {
            var name = e.slice(0, e.search("Scrib"));
            return $(name).toDataURL('image/png');
        });
    imgPNGs = entries;
}




function listRequest (e) {
    console.log("listRequest: " + e);
    var gbW = 800;
    var gbH = 700;
    var entries = gatherChecked('.item-box');

    listGBox(entries.slice(0,7).join(","), gbW, gbH);
}
