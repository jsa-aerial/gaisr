

var Templates = {
   marker: new Template (
       '<h2 style="margin-top: 0.75em;"><a href="#' + '#{name}"></a>' +
           '<input type="checkbox" value="on" name="#{name}" marker="#{taxid}" class="marker-box" onclick="toggleTaxon(this)">' +
           '<span class="marker">#{name}</span>' +
           '</h2>'),
   qentry: new Template (
      '<h2><input type="checkbox" value="on" name="#{name}" marker="#{taxon_id}" class="item-box">' +
           '<a href="#' + '#{name}">#{name}(#{sfcount})</a></h2>' +
      '<p>Name: #{name} (GID:#{gbid},V#{version})</p>' +
      '<p>Description: #{description}</p>' +
      '<p>Taxon: #{taxname} (ID:#{taxon_id})</p>' +
      '<p>#{ancestors}</p>' +
      '<p>#{sfcount} Listed Features</p>' +
      '<div class="links">' +
        '<h3>Links</h3>' +
        '<ul>' +
           '<li><a href="#" onclick="return createScriblGBWindow(\'#{name}\', \'#{name}\', #{bioentry_id}, #{gbid})">Scribl Features</a></li>' +
           '<li><a href="#" onclick="return createFeatureGBWindow(\'#{name}\', \'#{name}\', #{bioentry_id}, #{gbid})">List Features</a></li>' +
           '<li><a href="http://www.ncbi.nlm.nih.gov/nuccore/#{name}" onclick="return GB_showCenter(\'#{name}\', this.href, 700, 800)">NCBI</a></li>' +
           '<li><a href="http://www.ebi.ac.uk/genomes/bacteria.html" onclick="return GB_showCenter(\'#{name}\', this.href, 700, 800)">EMBL</a></li>' +
        '</ul>' +
      '</div>'),
   count: new Template (
      '#{count} Entries'),
   chat: new Template (
       '<p>Req: #{me}</p>' +
       '<p class="im-alice">SQL: #{alice}</p>'),

   ventry: new Template (
      '<h2><input type="checkbox" value="on" name="#{name}" marker="#{taxon_id}" class="item-box">' +
           '<a href="#' + '#{name}">#{name}(#{sfcount})</a></h2>' +
      '<p>Name: #{name} (ID:#{bioentry_id},V#{version})</p>' +
      '<p>Description: #{description}</p>' +
      '<p>Taxon: #{taxname} (ID:#{taxon_id})</p>' +
      '<p>#{ancestors}</p>' +
      '<p>Hit range: -#{delta}b / +#{delta}b</p>' +
      '<p>#{sfcount} Features within hit range</p>' +
      '<div class="links">' +
        '<h3>Links</h3>' +
        '<ul>' +
           '<li><a href="http://www.ncbi.nlm.nih.gov/nuccore/#{name}" onclick="return GB_showCenter(\'#{name}\', this.href, 700, 800)">NCBI</a></li>' +
           '<li><a href="http://www.ebi.ac.uk/genomes/bacteria.html" onclick="return GB_showCenter(\'#{name}\', this.href, 700, 800)">EMBL</a></li>' +
        '</ul>' +
      '</div>'),

   scribl: new Template (
      '<h3 class="feature-head">' +
         '<input type="checkbox" name="#{name}Scrib" class="scrib-item-box">' +
         '<span>#{name}[#{hitloc}]: EV #{evalue},</span>\
          <span id="#{name}fcount" style="display: none;">0</span>\
          <span class="scrib-box">Misc\
            <input id="#{name}mbx" name=#{name} type="checkbox" value="on" onclick="reScribl(this)"></span>\
          <span class="scrib-box">CDS\
            <input id="#{name}cbx" name=#{name} type="checkbox" value="on" onclick="reScribl(this)"></span>\
          <span class="scrib-box">Gene\
            <input id="#{name}gbx" name=#{name} type="checkbox" value="on" onclick="reScribl(this)"></span>\
       </h3>' +
      '<canvas id="#{name}" class="canvas"\
               width="#{scriblW}" height="#{scriblH}">\
       </canvas>'),
   scribcount: new Template (
       '#{feats} F, #{loci} L'),

   feature: new Template (
       '<h2><a href="#' + '#{sfid}">#{sfid}(#{sftype})</a></h2>' +
       '<div class="locs">' +
         '<h3>Locations</h3>' +
         '<ul id="#{sfid}loc"></ul>' +
       '</div>' +
       '<div class="nvs">' +
         '<h3>Qualifiers</h3>' +
         '<ul id="#{sfid}nvs"></ul>' +
       '</div>'),
   loc: new Template (
       '<p>#{start}..#{end}[#{len}], #{strand}</p>'),
   nv : new Template (
       '<p>#{name}: #{value}</p>'),

   reqHead: new Template (
       '<input type="text" name="act" value="#{act}" style="display: none;"/>' +
       '<input type="text" name="user" value="#{user}" style="display: none;"/>'),
   options1: new Template (
       '<div class="options" style="width: 21em; height: 29em;">' +
         '<fieldset id="options1">' +
           '<legend>Selections:</legend>' +
           '<input type="text" id="optSel" name="selections" style="display: none;"/>' +
           '<ul id="selUL"></ul>' +
         '</fieldset>' +
       '</div>'),
   options2: new Template (
       '<div class="options" style="width: 21em; height: 29em;">' +
         '<fieldset id="options2">' +
           '<legend>Options:</legend>' +
           '<ul id="optUL"></ul>' +
         '</fieldset>' +
       '</div>'),
   select: new Template (
       '<p>#{name}</p>')
}


function listGBox (title, w, h) {
    setupGBox(
        title, w, h,
        '<div id="Features" class="features" style="width: '+w+'px; height: '+h+'px;">\
         </div>',
        'Features',
        '<h3 class="feature-head">Features:\
          <span id="fcount" style="display: none;">0</span>\
          <img id="fspinner" src="ajax-aerial-small.gif" style="display: none;"/>\
         </h3>\
         <ul id="fresults" class="accordion">\
         </ul>');
}


function scriblGBox (title, w, h) {
    var scriblW = w - 25;
    var scriblH = h - 25;
    setupGBox(
      title, w, h,
      '<div id="Scribl" class="features"\
                        style="width: '+w+'px; height: '+h+'px;">\
       </div>',
      'Scribl',
      '<div id="scriblCtrl" class="features" style="display:none;">\
          <h3 class="feature-head">\
            <input type="checkbox" id="scriblCheckBox" value="on"\
                   onclick="scriblToggleChecks(this)"/>\
            Mark\
            <img id="scriblSpinner" src="ajax-aerial-small.gif"\
                 style="display: none;"/>\
            <div style="float: right;">\
              <input type="button" value="Names" onclick="genomeNamePage()"/> \
              <input type="button" value="PNG" onclick="genGenomePNG()"/> \
              <input type="button" value="SVG" style="display: none;"/>\
            </div>\
            <form id="scriblForm" class="actionForm"\
                  method="get" action="/mlab/action">\
              <input type="text" name="user" style="display: none;"/>\
              <input type="text" name="act" style="display: none;"/>\
              <input type="text" name="selections" style="display: none;"/>\
              <input type="text" name="filename" style="display: none;"/>\
            </form>\
          </h3>\
       </div>\
       <div id="genomeNamesDiv" style="display:none">\
         <textarea id="genomeNames" rows="20" style="width: 50%"></textarea>\
         <div style="display: inline; float: right;">\
           <form id="nameSaveForm" method="get" action="/mlab/action">\
             <fieldset>\
               <legend>DB Record Name</legend>\
               <input type="text" id="DBRec" name="DBRec"/>\
               <input type="button" value="Save2DB" onclick="dbSaveNames()"/>\
             </fieldset>\
           </form>\
           <form id="genFastaForm" method="get" action="/mlab/action">\
             <fieldset>\
               <legend>File Name</legend>\
               <input type="text" id="sfname" name="sfname"/>\
               <input type="button" value="Gen Fasta"\
                      onclick="scriblGenFasta()"/>\
             </fieldset>\
           </form>\
         </div>\
       </div>\
       <ul id="chartUL"></ul>');
}


