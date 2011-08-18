

function $jsa(element) {
  if (arguments.length > 1) {
    for (var i = 0, elements = [], length = arguments.length; i < length; i++)
      elements.push($(arguments[i]));
    return elements;
  }
  if (Object.isString(element))
    element = document.getElementById(element);
  return Element.extend(element);
}




new Ajax.Request('/aerial/action?p1=abc&p2=def',
                 {method: 'get',
                     onSuccess:
                   function(request) {
                     alert(request.responseText);
                     alert(request.responseXML);
                     alert(request.responseHeaders);
                   }
                 });


var poller = new Ajax.PeriodicalUpdater
  ('bucket', 'ajax.html',
   {method: 'get', insertion: 'bottom', fequency: 5, decay: 2});

poller.stop();
poller.start();


new Ajax.Updater({success: 'results', failure: 'error_log'},
		 '/aerial/qsite',
                 {method: 'get',
                  parameters: {food_type: 'waffles', taste: 'delicious'},
                  insertion: 'top'
                 });


new Ajax.Updater({success: 'results', failure: 'error_log'},
		 '/aerial/qsite',
                 {method: 'get',
                  parameters: {food_type: 'waffles'},
                  insertion: 'top'});

new Ajax.Updater({success: 'results', failure: 'error_log'},
		 '/aerial/qsite',
                 {method: 'get',
                  parameters: {food_type: 'waffles'},
                  insertion: 'top'});

var request = new Ajax.Request('/aerial/qsite',
			       {method: 'get',
				parameters:
				{json: 't', food_type: 'food', taste: 'xxx'},
				onSuccess: function (request) {
				 console.log(request.responseJSON);
				 }});

var tpl = new Template (
      '<h2><a href="#' + '#{id}">#{name}</a></h2>' +
      '<p>Title: #{title}</p>' +
      '<p>Dept: #{dept}</p>' +
      '<p>Company: #{company}</p>' +
      '<p>Office: #{location}</p>' +
      '<div class="links">' +
        '<h3>Profiles</h3>' +
        '<ul>' +
          '<li><a href="http://www.jigsaw.com">Jigsaw</a></li>' +
          '<li><a href="http://www.linkedin.com">LinkedIn</a></li>' +
        '</ul>' +
      '</div>');

var request = new Ajax.Request('/aerial/qsite',
			       {method: 'get',
				parameters:
				{json: 't', food_type: 'food', taste: 'xxx'},
				onSuccess: function (request) {
  alert(tpl.evaluate(request.responseJSON));
  console.log(request.responseJSON);
				 }});

var request = new Ajax.Request('/aerial/qsite',
			       {method: 'get',
				parameters:
				{json: 't', food_type: 'food', taste: 'xxx'},
				onSuccess: function (request) {
 var elt = $jsa(document.createElement('li'));
 elt.update(tpl.evaluate(request.responseJSON));
 alert(tpl.evaluate(request.responseJSON));
 console.log(request.responseJSON);}});



var request = new Ajax.Request('/aerial/qsite',
			       {method: 'get',
				parameters:
				{json: 't', food_type: 'food', taste: 'xxx'},
				onSuccess: function (request) {
 var elt = $jsa(document.createElement('li'));
 elt.update(Contacts._templates['contact'].evaluate(request.responseJSON));
 //alert(elt);
 $jsa('results').appendChild(elt);
 Accordion.init();
 //alert(tpl.evaluate(request.responseJSON));
 console.log(request.responseJSON);}});




$('results').down().next().next();

$('results').previous().up();

$('results').down('div'); == $('results').down().down().next();