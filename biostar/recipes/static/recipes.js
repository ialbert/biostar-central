/*
Code that handles recipe interface goes here.

Must call setup.js first load up functions used here.
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

// Submit edit request.
function submit_form() {


    var data = {
        'name': $('input[name=name]').val(),
        'uid': $('input[name=uid]').val(),
        'text': $('textarea[name=text]').val(),
        'image': $('input[name=image]').val(),
        'rank': $('input[name=rank]').val()
    };


    var url = '/recipe/edit/{0}/'.format(data.uid)

    $.ajax(url, {
            type: 'POST',
            dataType: 'json',
            data: data,
            success: function (data) {
                if (data.status === 'error') {
                    //popup_message($('#template'), data.msg, data.status);
                    alert(data.msg, data.status)
                    

                } else {
                    $('#message').html("OK")
                }

            },
            error: function (xhr, status, text) {

            }
        }
    );
}

// Toggles visible panels in recipe view.
function toggle_panels(elem_id, quick) {

    // Identify the element
    var elem = $(elem_id)

    // Move the element so it is first, thus always opens downwards.
    $(elem).parent().prepend(elem)

    //Hide all collapsible elements.
    $(".collapse").hide()

    // Open selected element.
    if (quick) {
        elem.show()
    } else {
        elem.show("slow", function () {
        });
    }

    // Set the window location hash
    window.location.hash = elem_id;

    // Remove active class on clickable object.
    $(".click").removeClass("active");

    // Formulate the selector.
    var selector = "[data-value='{0}']".format(elem_id)

    // Apply the active class to the selector.
    $(selector).addClass("active");

}

$(document).ready(function () {

    var script = prepare_codemirror($('#code textarea'), 700);
    var interface = prepare_codemirror($('#interface textarea'), 700);

    //script.refresh();
    //interface.refresh();

    // Select default open item.
    var open_id = window.location.hash || "#info";

    // Initial toggle has no animation.
    toggle_panels(open_id, 1)

    $(".click").click(function (event) {

        // Don't trigger other behaviors.
        event.preventDefault();

        // Find the targeted element is in the data-value of the clicked element.
        var target_id = $(this).data('value')

        // The selected page is already active.
        if (target_id === window.location.hash) {
            return;
        }

        // Toggle panels more slowly.
        toggle_panels(target_id, 0)

    });

    // Catch click on elements with submit types.
    $(":submit").click(function (event) {
        event.preventDefault();
        submit_form()
    });

});
