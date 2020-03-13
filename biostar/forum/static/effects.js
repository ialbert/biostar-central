
function allowDrop(ev) {
    ev.preventDefault();
    //alert(ev);
}

var children = '';
var dragged_over = '';

function drag(ev, elem_id) {
    //ev.preventDefault();
    ev.dataTransfer.setData("text", elem_id);
    //let elem = $('.droptarget[data-value="' + elem_id + '"]');

    //children = '{0}'.format(elem.data('thread'));

    ev.dataTransfer.setData("text", elem_id);
    ev.target.style.opacity = "0.4";

    //alert(elem_id);

}

function drag_leave(ev, elem) {

    //let elem = $('.droptarget[data-value="' + pid + '"]');
    elem.css('border', 'none');

    //elem.css('padding', '0');
    elem.css('opacity', '1');
    //elem.css('backgroundColor', 'none');
    dragged_over = '';
    //console.log(elem.data('value'), dragged_over);

}


function reset_drag(source, target) {
    console.log(target, source);
    source.css('opacity', '1');
    source.css('border', 'none');
    target.css('opacity', '1');
    target.css('border', 'none');

}

function drag_over(ev, elem) {
    dragged_over = '';
    //console.log(elem.attr('class'), elem.data('value'));
    dragged_over = elem.data('value');
    elem.css('border', '#c2ffc2 dashed 5px');

}


function drop(ev, elem_id) {

    ev.preventDefault();

    var source = ev.dataTransfer.getData("text");
    let source_elem = $('#' + source);
    source_elem.css('opacity', '1');
    let elem = $('.droptarget[data-value="' + dragged_over + '"]');

    //alert(dragged_over);
    $.ajax('/drag/and/drop/',
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'uid': source,
                'parent': dragged_over,

            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(elem, data.msg, data.status, 2000);
                    reset_drag(elem, source)
                } else {
                    source_elem.transition('zoom');
                    window.location.reload();
                    //window.location.href = data.redir;
                    popup_message(source_elem, "Moving Post", 'success', 1000);
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text);
                reset_drag(elem, source)
            }
        });

    ev.dataTransfer.clearData();
    dragged_over = '';

    reset_drag(elem, source)


}



function autocomplete_users(users_list) {
    // Add autocomplete to any text area element with autocomplete tag.
    var autocomplete = $('.autocomplete');

    // Map values in list to a list of dict [{key:'chosen', name:'displayed name'}....]
    var vals = $.map(users_list, function (value) {
        return {
            key: value,
            name: value
        };
    });

    function img_url(username) {
        let url = '/ajax/user/image/{0}/'.format(username);
        return "<li> <img class='ui circular image' style='display: inline' src={0} height='20' width='20' /> <b>{1}</b></li>".format(url, username);

    }

    // Autocomplete settings
    var AutocompleteSettings = {
        // Gets triggered at @
        at: "@",
        data: vals,
        displayTpl: img_url('${key}'),
        insertTpl: '@${key}',
        delay: 40
    };

    autocomplete.atwho(AutocompleteSettings);

}

$(document).ready(function () {

    $(this).on('ondragstart', '.draggable', function () {

    });

});

