/*
Code that handles recipe interface goes here.
*/

function prepare_codemirror(element, size) {

    var area = CodeMirror.fromTextArea(
        element[0],
        {
            lineNumbers: true,
            mode: 'shell',
        }
    );

    function update() {
        area.save();
    }

    area.on('change', update);

    area.setSize(null, size);

    return area
}

$(document).ready(function () {

    var script = prepare_codemirror($('#code textarea'), 700);
    var interface = prepare_codemirror($('#interface textarea'), 700);

    hash = window.location.hash || "#description" ;

    //script.refresh();
    //interface.refresh();

    // All collapsable elements
    collapse = $(".collapse")

    collapse.hide()

    //Show selected tab.
    $(hash).show()


    $(".item").click(function (event) {
        event.preventDefault();

        target = $(this)
        hash = target.data('value')

        console.log(target)
        console.log(hash)

        var current = window.location.hash;

        console.log(hash, current)

        if ('#' + hash === current){
            return;
        }

        window.location.hash = hash;

        collapse.hide("slow", function () {
            // Animation complete.
        });

        $("#" + hash ).show("slow", function () {
            // Animation complete.
            //window.location.hash = hash;
        });

    });

});
