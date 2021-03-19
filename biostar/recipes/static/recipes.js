/*
Code that handles recipe interface.

Must call setup.js first initialize additional functionality used here.

*/


function codemirror_callbacks() {

    return {
        'Shift-Enter': (cm) => {
            update_preview();
            popup_message($("#interface"), "Updated the interface preview", "success", 4000)
        }
    }
}

function init_codemirror(element, size) {

    var area = CodeMirror.fromTextArea(
        element[0],
        {
            autoRefresh: true,
            lineNumbers: true,
            mode: 'shell',
            theme: 'idea',
            extraKeys: codemirror_callbacks(),

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
function remove_messages() {
    $('.notify').remove()
}


function preview_template(fields) {

    return '<div class="ui segment run"><form class="ui form">'+
        '<span class="preview">' + fields+ '</span>'+
        '<div class="field">' +
        '<button type="submit" class="ui green disabled button">' +
        '    <i class="check icon"></i>Run' +
        '</button>' +
        '<a class="ui disabled button">' +
        '    <i class="redo icon"></i>Cancel' +
        '</a>' +
        '</div></form></div>'
}


function update_preview() {

    let toml = $('#interface_editor').val();
    let project = $('#interface').closest('.grid').data("project");
    let url = '/preview/json/';
    var id = get_id();

    $.ajax(url, {
        type: 'POST',
        dataType: 'json',
        data: {
            'recipe': id,
            'toml': toml
        },
        success: function (data) {
            if (data.status === 'error') {
                popup_message($("#interface"), data.msg, data.status, 5000);
                return
            }
            $('#preview').html(preview_template(data.html));
            //alert(preview_template(data.html));
            //$('#preview').show()
        }
    });
}

// Submits a form request.
function submit_form(elem) {

    // The form the button belongs to.
    form = elem.closest("form");

    // Get the recipe id.
    var id = get_id()

    // Bind the closest form data.
    var data = new FormData(form.get(0));

    // Recipe id must be used here.
    var url = '/ajax/recipe/edit/{0}/'.format(id)

    // Remove any prior notification that may exist.
    remove_messages();

    $.ajax(url, {
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            data: data,
            success: function (resp) {
                if (resp.status === 'error') {
                    show_message(elem, resp.msg, "error");
                    //flash("fadeout_error")
                } else {
                    if (window.location.hash === '#edit') {
                        toggle_panels('#info', 1);
                    } else {
                        //flash("fadeout_final")

                        popover_message(form, resp.msg, "success")
                    }
                    update_panels();
                }
            },
            error: function (xhr, status, error) {
                var text = "Error: status={0} error={1}".format(status, error)
                show_message(elem, text, "error")
            }
        }
    );
}

function flash(cls) {
    elem = $(".CodeMirror");

    function fadeout() {
        elem.removeClass(cls);
        elem.addClass("fadeout_final")
    }

    elem.removeClass("fadeout_final")
    elem.addClass(cls);

    setTimeout(fadeout, 1500);

}

// Toggles visible panels in recipe view.
function toggle_panels(elem_id, quick) {

    // Identify the element
    var elem = $(elem_id);

    // Move the element so it is first, thus always opens downwards.
    $(elem).prepend(elem);

    //Hide all collapsible elements.
    $(".collapse").hide();

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
    var selector = "[data-value='{0}']".format(elem_id);

    // Apply the active class to the selector.
    $(selector).addClass("active");


}

// Updates content in dynamic panels
function update_panels() {

    var panels = ['info', 'run', 'results', 'details'];
    var server = "/get/part/{0}/{1}/";
    var id = get_id();

    function loader(name) {
        var node = $('#{0}'.format(name));
        var url = server.format(name, id);
        node.load(url, function (response, status, xhr) {

        });
    }

    panels.forEach(loader);
}



function validate_toml(input_field, toml_text){
    /*
     Check if a new toml field is already present in the text.
     Increments the field key to resolve conflicts.
     */

    // List of the current toml parameter being added
    var input_list = input_field.trim().split("\n");

    /// Get the key of this interface parameter
    var current = input_list[0].trim();
    var key = current.slice(1, -1);
    var i = 1;

    // Loop until we get to a unique key
    while (toml_text.includes(current, 0)){
        new_key = "{0}_{1}".format(key, i);
        current = "[{0}]".format(new_key);
        console.log(current);
        i ++;
    }

    input_list[0] = current;

    // Join the new toml
    let input_str = "\n\n" + input_list.join("\n");

    return input_str


}


// Binds events dynamically.
function bind_events() {

    $(document).on("click", '.click', function (event) {
        event.preventDefault();

        // Find the targeted element is in the data-value of the clicked element.
        var target_id = $(this).data('value');

        // The selected page is already active.
        if (target_id === window.location.hash) {
            return;
        }

        // Toggle panels slowly.
        toggle_panels(target_id, 0)
    });

    // Forms with ajax submissions.
    $(":submit.ajax").click(function (event) {
        event.preventDefault();
        submit_form($(this));
        // Update the preview on form submit.
        update_preview();

    });

}

// The recipe id obtained from the page.
function get_id() {
    return $("#recipe_id")[0].value
}

$(document).ready(function () {

    // Select default open item.
    var open_id = window.location.hash;

    // Update information in dynamic panels when id is provided.
    if (open_id){
        // Initial toggle has no animation.
        toggle_panels(open_id, 1);
    }else{
        // Show the info panel by default when a hash is not set.
        $(".collapse:not(#info)").hide();
    }

    update_panels();

    // Bind the events.
    bind_events();

});
