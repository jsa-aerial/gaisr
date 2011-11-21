

var Accordion = {
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
                        var kthChild = folds[j].childNodes[k];
                        if (kthChild.nodeName.toLowerCase() != "h2") {
                            accordionContent.appendChild(kthChild);
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
                var activeFold = $(location.hash.substring(1));
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



