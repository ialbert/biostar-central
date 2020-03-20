
function drop(ev, elem) {

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

function bind_drag_and_drop() {
    //var dragging = '';
    $(".post > .body > .content > .droppable").droppable(
        {
            accept: ".post",
            drop: function (event, ui) {
                //alert($(this).closest(".post").data("value"));
                //alert(ui.draggable.data("value"));
                //drop(event, $(this));
            },

        });

    $('.draggable').mousedown(function () {
        $(this).css('cursor', 'grabbing');
        var post = $(this).closest('.post');
        post.draggable(
        {
            addClasses: false,
            scroll: false,
            helper: 'clone',
            iframeFix: true,
            opacity: .7,
            containment: $('body'),
            revert:true,
            zIndex: 100,
            cursor:'grabbing'

        });
        //post.draggable( "disable" );
    });



}
