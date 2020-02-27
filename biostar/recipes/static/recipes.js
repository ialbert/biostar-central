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

// Toggles visible panels in recipe view.
function toggle_panels(elem_id, quick) {

    // Identify the element
    var elem = $(elem_id)

    // Move the element so it is first, thus always opens downwards.
    $(elem).parent().prepend(elem)

    // Select collapsible elements
    var collapse = $(".collapse")

    //Hide should always be fast
    collapse.hide()

    // Close all collapsible elements.
    if (quick) {
        elem.show()
    } else {
        elem.show("slow", function () {
            // Animation complete.
        });
    }

    // Set the window location hash
    window.location.hash = elem_id;

    // Remove active class on clickable object.
    $(".click").removeClass("active")

    // Highlight the selector.
    $('[data-value="' + elem_id + '"]').addClass("active")

}

$(document).ready(function () {

    var script = prepare_codemirror($('#code textarea'), 700);
    var interface = prepare_codemirror($('#interface textarea'), 700);

    //script.refresh();
    //interface.refresh();

    // Select default open item.
    var open_id = window.location.hash || "#description";

    // Initial toggle has no animation.
    toggle_panels(open_id, 1)

    $(".click").click(function (event) {

        // Don't trigger other behaviors.
        event.preventDefault();

        // The clicked element
        var elem = $(this)

        // Find the targeted element.
        var target_id = elem.data('value')

        // The selected page is already active.
        if (target_id === window.location.hash) {
            return;
        }

        // Toggle panels more slowly.
        toggle_panels(target_id, 0)

    });

});
