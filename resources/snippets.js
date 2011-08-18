/*

  GB_showCenter('testing...', '/meyer-lab.html', 650, 750); // Open a window
  gbw = $('GB_window');
  GBiframe = $$('.GB_frame')[0]; // More direct...
  GBiframe = gbw.children[1].children[0];
  GB_hide(); // NOTE this is global!

  ;; Snippets of various greybox stuff.

  ;; This one supposedly closes the window programatically (here via a
  timer but could simply be a direct call.)
  
  <body  onload="self.setTimeout('parent.parent.GB_hide()', 5000);">

  ;; This one supposedly supports calling a form from within greybox:

  onclick="return GB_show('Google', this.href)"

  so adding this to your button or input should work..
  like this:

  <input type="submit" value="submit" onclick="return GB_show('Google', this.href)" />

  However, with that submit button you don't have an href attribute to
  get the url.  Also I would recommend a button so that the page
  doesnt try to postback.

  <input type="button/or submit" value="Some title/label"
         onclick="return GB_showCenter('Page title', 'http://www.google.com')"/>



  function clickButton(url){
  // create the queryString out of the form --> replace FORMID with the
  // ID of the form you want to submit but be sure to leave the #
  var queryString = "?" + $("#FORMID").serialize();

  // pass your form data querystring to your server side code and pop the
  // greybox with the data that is returned
  return GB_showCenter("Greybox Title", url + queryString);
}
</script>

And your link will look similar to this:
<a href="#" onclick="clickButton('http://youraddress.com/serverPage.asp')">Click Here</a>



*/


function showDialog(t, reload){
    if(typeof t != 'string') {
        t = $(this).title() || $(this).text() || $(this).href();
    }
    var callback = function(){
        $('form', $('#GB_frame').get()).ajaxForm(
            {'target': '#GB_frame',
             'after': function() {
                 var status = $('#status');
                 // yeah, I use a class name to signal if I'm done. YGAPWT?
                 var code = status.attr('class');
                 var status_msg = status.html()
                 var content = $('#GB_frame').html(); // the whole page in the dialog
                 if (code != 'ok') {
                      // call myself again, with the content of the dialog box
                     showDialog(t, content);
                 }
                 else {
                     $.GB_hide();
                     alert(status_msg);
                 }
             } } );
    }; // end of callback.
    if (!reload) { // we were called for the first time; create a dialog box!
        var url = $(this).href();
        var arguments = null;
        $.GB_show('about:blank', {
            height: 400,
            width: 400,
            animation: true,
            overlay_clickable: true,
            caption: t
        });
        // We don't want the iframe greybox gives us:
        $('#GB_frame').remove();
        $("#GB_window").append("<div id='GB_frame'></div>");
        $("#GB_frame").load(url, // URL
                            arguments, // Params
                            callback );
    }
    else { // we were called again, because we're not done yet: just reload the HTML
        $("#GB_frame").html(reload);
        callback();
    }
    return false;
}

$(document).ready(function(){ $("a.dialog").click(showDialog); });

