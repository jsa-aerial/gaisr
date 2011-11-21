

/*
  --- Action Request Processing ---
*/


// Histories are just typical circular buffers over an array.  Since
// JS arrays are expandable, we don't need to allocate full size up
// front.
var inputHist = {
    val: [],
    index: 0,
    len: 0,
    max: 50,
    lastKey: undefined
}

function initHistory () {
    var q = Cookies.qhist;
    inputHist.val = (q) ? q.split(',') : [];
    inputHist.len = inputHist.val.length;
    inputHist.index = 0;
    //console.log(inputHist.val);
}


function handleInputHistory (e) {
    var element = Event.element(e); // element = input "query" text box
    var key = e.keyCode;
    var lastKey = inputHist.lastKey;
    var len = inputHist.len;

    // This looks "upside down", because push and pop of JS arrays
    // adds to the end of the array and there is no list type.
    switch (key) {
        case Event.KEY_UP:
            e.stop;
            if (lastKey == Event.KEY_DOWN) {
                inputHist.index = (inputHist.index - 1 + len) % len;
            }
            inputHist.index = (inputHist.index - 1 + len) % len;
            element.value = inputHist.val[inputHist.index];
            break;

        case Event.KEY_DOWN:
            e.stop;
            if (lastKey == Event.KEY_UP) {
                inputHist.index = (inputHist.index + 1) % len;
            }
            element.value = inputHist.val[inputHist.index];
            inputHist.index = (inputHist.index + 1) % len;
            break;
    }

    inputHist.lastKey = key;
}


// Update history and reset index
function addHistoryItem (item) {
    if (inputHist.len < inputHist.max) {
        inputHist.val.push(item);
        inputHist.len = inputHist.val.length;
    } else {
        inputHist.val.shift(1);   // Drop oldest
        inputHist.val.push(item); // len stays max
    }
    inputHist.index = 0;
    setCookie("qhist", inputHist.val.toString(), 365);
    return item;
}




function evaluePass (e) {
    var cutoff = $('evinput').getValue();
    //console.log(hitFeatures[e]["evalue"], " < ", cutoff,
    //            hitFeatures[e]["evalue"] < cutoff);
    return hitFeatures[e]["evalue"] < cutoff;
}


// Global selection and deselection of all current query results
//
function toggleChecks (elt) {
    console.log("Results checked -- " + elt.getValue());
    var boxes = $$('.item-box');
    var taxonBoxes = $$('.marker-box');
    var evcheck = $('evspan').visible();
    var chkcount = 0;
    if (elt.getValue() == null) {
        boxes.each(function (b) {b.setValue(false);});
        taxonBoxes.each(function(b){b.setValue(false);});
    } else {
        boxes.each(function (b) {
                if (evcheck) {
                    if (evaluePass(b.name)) {
                        b.setValue(true);
                        chkcount = chkcount + 1;
                    }
                } else {
                    b.setValue(true);
                }});
    }
    if (evcheck) {
        $('evcutoffcnt').update(chkcount);
    }
}


// Selection on taxon category
//
function toggleTaxon (elt) {
    console.log("Taxon " + elt.getAttribute('name') +
                " checked -- " + elt.getValue());
    var taxid = elt.getAttribute('marker');
    var boxes = $$('.item-box').findAll(function(x){
        return x.getAttribute('marker') == taxid;});
    if (elt.getValue() == null) {
        boxes.each(function(b) {b.setValue(false);});
    } else {
        boxes.each(function(b) {b.setValue(true);});
    }
}


// Gather up the names/ids of all the selected (i.e., checked) items
// in the current results set.  This is used in the actonSelected
// primary function for handling action requests.
//
function gatherChecked (boxClass) {
  return $$(boxClass).filter(function (e) {
      return e.getValue() == "on";
  }).pluck('name');
}



// OK, it is not clear that submitting the content of the action form
// with secondary events makes all that much sense.  It may make
// rather more sense to simply call the observer function (see
// bindForm in meyer-lab.js) directly (after first factoring it out as
// a named function).  Events decouple things more, but it isn't clear
// that here that is of any real use.  The downside is that, well, it
// decouples things more.
//
// The action form contains a text field called SELECTIONS, and a text
// field called ACTION, which are what actually gets submitted as the
// serialized parameter(s).  The form also contains a number of
// buttons which act as "submit drivers".  This lets us use only one
// form and only the two input fields for the action and selections
// (instead of one set per button).
//
function fireRequestForm (action) {
    $('actionForm').fire('action:submit', {action: action});
}


// Following set make action specific to that requested
//
function jobRequest (e) {
    console.log("Job Request");
    fireRequestForm("job");
}

function blastRequest (e) {
    console.log("Blast Request");
    fireRequestForm("blast");
}

function cmfindRequest (e) {
    console.log("CMFind Request");
    fireRequestForm("cmfind");
}

function dbBuildRequest (e) {
    console.log("DBbuild Request");
    fireRequestForm("dbbuild");
}

function genFastaRequest (e) {
    console.log("genFastaRequest...");
    fireRequestForm("genfasta");
}

function cmsearchRequest (e) {
    console.log("CMSearch Request");
    fireRequestForm("cmsearch");
}




var actionInfo = {
    reqCallBackMap: {blast: processBlastRequest,
                     cmfind: processCMFindRequest,
                     dbbuild: processDBbuildRequest,
                     genfasta: processGenFastaRequest,
                     cmsearch: processCMSearchRequest},
    blastOptions: {
        formFields: ['<span class="options">' +
                     'Type:' +
                     '<input type="radio" name="pgm" value="1"\
                             class="suboptions" checked>blastn' +
                     '<input type="radio" name="pgm" value="2"\
                             class="suboptions">tblastn' +
                     '</span>',

                     '<span class="options">' +
                     'Strand:' +
                     '<input type="radio" name="strand" value="1"\
                             class="suboptions" checked>plus' +
                     '<input type="radio" name="strand" value="2"\
                             class="suboptions">minus' +
                     '</span>',

                     '<span class="options">' +
                     'Evalue:' +
                     '<input type="text" name="evalue" size="1" value="10"\
                             class="suboptions">' +
                     '</span>',

                     '<span class="options">' +
                     'Word Size:' +
                     '<input type="text" name="wordsize" size="1" value="8"\
                             class="suboptions">' +
                     '</span>']
    },

    resCallBackMap: {blast: processBlastResults,
                     cmfind: processCMFindResults,
                     dbbuild: processDBbuildResults,
                     genfasta: processGenFastaResults,
                     cmsearch: processCMSearchResults}
}


function processBlastResults(request) {
    console.log("processBlastResults");
}

function processCMFindResults(request) {
    console.log("processCMFindResults");
}

function processDBbuildResults(request) {
    console.log("processDBbuildResults");
}

function processGenFastaResults(request) {
    console.log("processGenFastaResults");
}

function processCMSearchResults(request) {
    console.log("processCMSearchResults");
}


function submitAction (e) {
    $('actionOptions').fire('action:do');
}

function doAction (form, e) {
    console.log("doAction");
    var opt1val = $('selUL');
    var s = gatherChecked('.item-box');
    //var s = opt1val.childElements().map
    //    (function (e) {return e.down('p').textContent;});
    var optSel = $('optSel');
    optSel.setValue(s.join(" "));
    console.log(form.serialize());
}




function processBlastRequest(request) {
    console.log("processBlastResults");
    var opt2 = $('optUL');
    var fields = actionInfo.blastOptions.formFields;
    fields.each(function(f){
            var elt = $(document.createElement('li'));
            elt.update(f);
            opt2.appendChild(elt);
        });
}

function processCMFindRequest(request) {
    console.log("processCMFindRequest");
}

function processDBbuildRequest(request) {
    console.log("processDBbuildRequest");
}

function processGenFastaRequest(request) {
    console.log("processGenFastaRequest");
}

function processCMSearchRequest(request) {
    console.log("processCMSearchRequest");
}


function processActionRequest (request) {
    var items = request.responseJSON;
    var headMap = items[0];
    var act = headMap.action.toLowerCase();
    var user = $('user').getValue();
    var reqHead = {act: act, user: user};
    var actOptions = $('actionOptions');

    actOptions.insert(Templates['reqHead'].evaluate(reqHead));
    actOptions.insert(Templates['options1'].evaluate());
    actOptions.insert(Templates['options2'].evaluate());
    actionInfo.reqCallBackMap[act](request);

    var submitButton = $('optionSubmit');
    submitButton.observe('click', submitAction);

    actOptions.observe('action:do', function (e) {
      e.stop();
      doAction(this, e);
    });

    var opt_ul = $('selUL');
    items.each(function (item, index) {
        if (index > 0) {
            var elt = $(document.createElement('li'));
            elt.update(Templates['select'].evaluate(item));
            opt_ul.appendChild(elt);
        }
    });
}


function actonSelected(form, e) {
    var s = gatherChecked('.item-box');
    var a = e.memo.action;
    var u = $('user').getValue();

    if (s.length == 0) {
        alert(a + ' requires selections');
        return undefined;
    }

    if (u == "") {
        u = prompt('Action ' + a + ' missing user', "enter user name here");
        if (u == "") {
            alert(a + 'requires user');
            return undefined;
        }
        setUser(u);
    }

    var user = form.children[0];
    var act = form.children[1];
    var txt = form.children[2];
    user.setValue(u);
    act.setValue(a);
    txt.setValue(s.join(" "))

    // Setup GB window for request...
    //
    setupGBox(
        a, 600, 700,
        '<div id="Request" class="action" style="width: 600px; height: 700px;">\
         </div>',
        'Request',
        '<form id="actionOptions" method="post" action="/mlab/doaction">\
         <h3 class="action-head">' + a.capitalize() + ':' +
        '<input type="button" id="optionSubmit" value="Submit"\
                style="margin-left: 1em;"/>' +
        '</h3>' +
        '</form>'
    );

    var launchAjax = function () {
        new Ajax.Request(form.action, {
            method: 'post',
            parameters: form.serialize(),
            onSuccess: function (request) {
                processLog(request);
                processActionRequest(request);
            }
        });
    }
    // Faked sleep - gives the greybox time to fully paint before Ajax
    // response takes over
    setTimeout(launchAjax, 500);
}
