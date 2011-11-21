

// "Main" stuff...

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
    console.log('GAISR client init...');
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
    console.log('GAISR client initialized.');
  });

