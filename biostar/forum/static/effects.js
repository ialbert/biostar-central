
function allowDrop(ev) {
    ev.preventDefault();
}

var dragged_over = '';


function drag(ev, elem) {
    var post = elem.closest('.post');
    var uid = post.data('value');
    ev.dataTransfer.setData("text", uid);
    post.css('opacity', "0.4");

}

function drag_leave(elem) {
    var post = elem.closest('.post');
    post.css('border', 'none');
    post.css('opacity', '1');
    dragged_over = '';

}


function reset_drag(source, target) {
    console.log(target, source);
    source.closest('.post').css('opacity', '1');
    source.closest('.post').css('border', 'none');
    target.closest('.post').css('opacity', '1');
    target.closest('.post').css('border', 'none');

}

function drag_over(elem) {
    dragged_over = '';
    var post = elem.closest('.post');
    //console.log(elem.attr('class'), elem.data('value'));
    dragged_over = post.data('value');
    post.css('border', '#c2ffc2 dashed 5px');

}


function drop(ev, elem) {

    ev.preventDefault();

    var source = ev.dataTransfer.getData("text");
    let source_elem = $('#' + source);
    source_elem.css('opacity', '1');
    //let elem = $('.droptarget[data-value="' + dragged_over + '"]');

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

    $('.draggable').each(function () {
       $(this).closest('.post').attr('draggable', true);
    });
    $(this).on('ondragover', '.droppable', function (event) {
        allowDrop(event)
    });
    $(this).on('ondrop', '.droppable', function (event) {
        drop(event, $(this));
    });
    $(this).on('ondragstart', '.draggable', function (event) {
        drag(event, $(this))
    });

    $(this).on('ondragover', '.draggable', function (event) {
        drag_over($(this))
    });

    $(this).on('ondragleave', '.draggable', function (event) {
        drag_leave($(this))
    });

});

