window.onload=function() {
    var preComEL= document.getElementById('pre-compiled-radio');
    var selfComEL = document.getElementById('self-compiled-radio');

    preComEL.addEventListener('change', function() {
        toggleField("pre-compiled-more-info", "self-compiled-more-info");
    });

    selfComEL.addEventListener('change', function() {
        toggleField("self-compiled-more-info", "pre-compiled-more-info");
    });
}

function toggleField(id1, id2) {
    var elem1 = document.getElementById(id1);
    var elem2 = document.getElementById(id2);
    elem1.removeAttribute("hidden");
    elem2.setAttribute("hidden", null);
};

