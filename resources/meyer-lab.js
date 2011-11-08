

var Accordion =
{
  init: function(accordList) {
    Accordion.frameRate = 25;
    Accordion.duration = 0.5;

    var accordions = accordList || $$(".accordion");

    for (var i = 0; i < accordions.length; i++) {
      var folds = accordions[i].childNodes;

      for (var j = 0; j < folds.length; j++) {

        if (folds[j].nodeType == 1) {
          var accordionContent = document.createElement("div");
          accordionContent.className = "accordionContent";

          for (var k = 0; k < folds[j].childNodes.length; k++) {
            if (folds[j].childNodes[k].nodeName.toLowerCase() != "h2") {
              accordionContent.appendChild(folds[j].childNodes[k]);
              k--;
            }
          }

          folds[j].appendChild(accordionContent);
          folds[j]._accordionContent = accordionContent;

          Accordion.collapse(folds[j]);
          var foldLinks = folds[j].getElementsByTagName("a");
          var foldTitleLink = foldLinks[0];
          foldTitleLink.observe('click', Accordion.clickListener);

          for (var k = 1; k < foldLinks.length; k++) {
            foldLinks[k].observe('focus', Accordion.focusListener);
          }
        }
      }

      if (location.hash.length > 1) {
        var activeFold = document.getElementById(location.hash.substring(1));
        if (activeFold && activeFold.parentNode == accordions[i]) {
          Accordion.expand(activeFold);
        }
      }
    }
  },


  collapse: function(fold) {
    var content = fold._accordionContent;
    content._height = parseInt(content.style.height, 10);
    content._increment = content._height / (Accordion.frameRate *
                                            Accordion.duration);
    if (fold.hasClassName('expanded')) {
      clearTimeout(content._timer);
      Accordion.collapseAnimate(content);
    } else {
      fold.addClassName('collapsed');
    }
  },


  collapseAnimate: function(content) {
    var newHeight = content._height - content._increment;

    if (newHeight < 0) {
      newHeight = 0;
      content.parentNode.removeClassName('expanded');
      content.parentNode.addClassName('collapsed');
    } else {
      content._timer = setTimeout(function() {
          Accordion.collapseAnimate(content);
        }, 1000 / Accordion.frameRate);
    }

    content._height = newHeight;
    content.style.height = Math.round(newHeight) + "px";
  },


  collapseAll: function(accordion) {
    var folds = accordion.childNodes;
    for (var i = 0; i < folds.length; i++) {
      if (folds[i].nodeType == 1) {
        Accordion.collapse(folds[i]);
      }
    }
  },


  expand: function(fold) {
    var content = fold._accordionContent;
    Accordion.collapseAll(fold.parentNode);

    if (!fold.hasClassName('expanded')) {
      content.style.height = "0";
      content._height = 0;
      fold.removeClassName('collapsed');
      fold.addClassName('expanded');
      content._increment = content.scrollHeight / (Accordion.frameRate *
                                                   Accordion.duration);
      Accordion.expandAnimate(content);
    }
  },


  expandAnimate: function(content) {
    var newHeight = content._height + content._increment;
    if (newHeight > content.scrollHeight) {
      newHeight = content.scrollHeight;
    } else {
      content._timer = setTimeout(function() {
          Accordion.expandAnimate(content);
        }, 1000 / Accordion.frameRate);
    }

    content._height = newHeight;
    content.style.height = Math.round(newHeight) + "px";
    content.scrollTop = 0;
  },


  clickListener: function(event) {
    var fold = this.parentNode.parentNode;
    if (fold.hasClassName('collapsed')) {
      Accordion.expand(fold);
    } else {
      Accordion.collapse(fold);
    }
    event.preventDefault();
  },


  focusListener: function(event) {
    var element = this;
    while (element.parentNode) {
      if (element.parentNode.className == "accordion") {
        Accordion.expand(element);
        return;
      }
      element = element.parentNode;
    }
  }
};




var EXX = "";

function bindForm() {
  $('searchForm').observe('submit', function(e) {
    // Stop the browser from trying to handle this
    e.stop();
    // Send off AJAX request to server.
    // Set processResults as watch on result event.
    makeQuery(this, e);
  });

  $('query').observe('keyup', handleInputHistory);

  var actionForm = $('actionForm');

  actionForm.children[3].observe('click', jobRequest);
  actionForm.children[4].observe('click', blastRequest);
  actionForm.children[5].observe('click', cmfindRequest);
  actionForm.children[6].observe('click', dbBuildRequest);
  actionForm.children[7].observe('click', cmsearchRequest);

  actionForm.observe('action:submit', function(e) {
    // Stop the browser from trying to handle this event,
    // actionSelected will gather selections and send Ajax request...
    e.stop();
    //console.log("E = " + e + " Action: " + e.memo.action);
    //EXX = e;
    actonSelected(this, e);
  });

  var viewDiv = $('view-div');
  viewDiv.children[0].observe('click', scriblRequest);
  viewDiv.children[1].observe('click', listRequest);
}




document.observe("dom:loaded", function () {
    var u = Cookies.user;
    if (u) $('user').setValue(u);
    initHistory();

    bindForm();
    Accordion.init($$('#results'));
    document.observe("results:updated", populateResults);

    var navList = $('navList');
    navList.childElements().slice(1).each
        (function (x) {x.down('a').addClassName('tab-inactive');});
    navList.children[0].down('a').addClassName('tab-active');
    activeTab = navList.children[0].down('a');

    Ajax.Responders.register({
      onCreate: function () {
          if ($('fspinner')) {
            $('fcount').hide();
            $('fspinner').show();
          } else {
            $('count').hide();
            $('spinner').show();
          }
      },
      onComplete: function () {
          console.log("***Ajax COMPLETED");
          if (0 == Ajax.activeRequestCount) {
            if ($('fspinner')) {
                $('fspinner').hide();
                $('fcount').show();
            } else {
                $('spinner').hide();
                $('count').show();
            }
          };
        }
      });
  });



