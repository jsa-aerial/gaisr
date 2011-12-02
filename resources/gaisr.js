
// A handful of general "utility" functions up front that are used to
// include files and used in included files.  No real good place to
// put these...


// Used to include other js files/libs.  Unclear on whether this is
// really sufficient, but seems to work for our use...
//
function includeJS(filenames) {
    filenames = (typeof(filenames) != "object") ? [ filenames ] : filenames;
    if (filenames.length > 0) {
        var file = filenames[0];
        var head = document.getElementsByTagName('head')[0];
        var elt = document.createElement('script');
        elt.src = file;
        elt.type = 'text/javascript';

        var nextfn = function() {
            includeJS(filenames.slice(1));
        }
        elt.onload = nextfn;
        //console.log("Load '" + file + "', " + filenames);
        head.appendChild(elt);
    }
}


// Under some (totally uncharacterized circumstances) prototype $ may
// not work _at the FireBug REPL_.  Some(how/one) monkey patches over
// it and it's not us!!
//
function gaisr$(element) {
  if (arguments.length > 1) {
    for (var i = 0, elements = [], length = arguments.length; i < length; i++)
      elements.push($(arguments[i]));
    return elements;
  }
  if (Object.isString(element))
    element = document.getElementById(element);
  return Element.extend(element);
}


// Basically a simple "logging" debug "wrapper" for Ajax response
//
var reqJson = "NA";

function processLog(request) {
    console.log("START FOO")
    reqJson = request;
    console.log("DONE FOO");
}


function stopEvent (e) {
    e.stop();
    e.cancelBubble = true;
    if (e.stopPropagation) e.stopPropagation();
}


// These are redefined later, but needed to placate Chrome (and others?)
function uploadDone () {}
DOMLoaded = undefined;

// Pull in the rest of the client...
//
var filenames = ['JSLibs/prototype.js',
                 'JSLibs/scriptaculous.js',
                 'JSLibs/Scribl/Scribl.min.js',
                 'JSLibs/Scribl/Scribl.svg.js',
                 'cookies.js',
                 'html-templates.js',
                 'entry-info.js',
                 'gaisr-actions.js',
                 'gaisr-effects.js',
                 'gaisr-main.js'];

includeJS(filenames);


var gaisrInit = false;
var gaisrWin = undefined;
function syncDOMLoad (ms) {
    if (DOMLoaded == undefined) {
        setTimeout(syncDOMLoad, ms);
    } else if (gaisrInit == false) {
	//	if (gaisrWin != window) {
	//    console.log("gaisrWin: " + gaisrWin + ", /= " + window);
	//    gaisrWin = window;
	//}
	//console.log("Bef: gaisrInit = '" + gaisrInit + "'");
        gaisrInit = true;
        DOMLoaded();
	//console.log("Aft: gaisrInit = '" + gaisrInit + "'");
    } else {
        console.log("onload fired twice - init once");
    }
}

onload = function () {syncDOMLoad(500);}

