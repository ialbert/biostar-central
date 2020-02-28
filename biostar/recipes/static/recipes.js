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

// Removes notifications
function remove_notify() {
    $('.notify').remove()
}

// Shows messages.
function show_message(elem, mesg, status) {
    remove_notify()
    var node = $("<div class='notify'>{0}</div>".format(mesg));
    node.addClass("ui {0} message popover".format(status));
    elem.closest("form").addClass(status).prepend(node)
    node.delay(1000).fadeOut(500, function () {
        $(this).remove()
    });
}

function popup_message(elem, message, cls, timeout) {
    timeout = typeof timeout !== 'undefined' ? timeout : 1000;

    elem = elem.find("textarea")

    var text = $('<div class="popover"></div>');
    var tag = $(text).insertBefore(elem)
    tag.addClass(cls)
    tag.text(message)
    tag.delay(timeout).fadeOut(500, function () {
        $(this).remove()
    });
}

function flash(cls){
    elem = $(".CodeMirror")
    function explode(){
        elem.removeClass(cls)
        elem.addClass("fadeout")
    }
    setTimeout(explode, 1000);
    elem.removeClass("fadeout")
    elem.addClass(cls)
}


// Submits a form request.
function submit_form(elem) {

    // The form the button belongs to.
    form = elem.closest("form");

    // Get the data from the form.
    var data = {
        'json_text': form.find("textarea[name=json_text]").val() || '',
        'template': form.find("textarea[name=template]").val() || '',
        // This variable is special and is used as submit id.
        'id': form.find('input[name=id]').val()
    };

    // Recipe id must be used here.
    var url = '/recipe/ajax/edit/{0}/'.format(data.id)
    
    $.ajax(url, {
            type: 'POST',
            dataType: 'json',
            //processData: false,
            data: data,
            success: function (res) {
                if (res.status === 'error') {
                    show_message(elem, res.msg, res.status)
                } else {
                    //show_message(elem, res.msg, "success")
                    flash("fadein_success")
                }
            },
            error: function (xhr, status, error) {
                var text = "Ajax error: status={0} error={1}".format(status, error)
                show_message(elem, text, "error")
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
    var interface = prepare_codemirror($('#interface textarea'), 400);

    //script.refresh();
    //interface.refresh();

    // Select default open item.
    var open_id = window.location.hash || "#info";

    // Initial toggle has no animation.
    toggle_panels(open_id, 1);

    $(".click").click(function (event) {

        // Don't trigger other behaviors.
        event.preventDefault();

        // Find the targeted element is in the data-value of the clicked element.
        var target_id = $(this).data('value');

        // The selected page is already active.
        if (target_id === window.location.hash) {
            return;
        }

        // Toggle panels more slowly.
        toggle_panels(target_id, 0)

    });

    // Catch click on elements with submit types.
    $(":submit.ajax").click(function (event) {
        event.preventDefault();
        submit_form($(this))
    });

});
