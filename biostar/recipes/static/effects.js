


function move_object(parent_elem, source_elem, child_elem){

    var source_id = source_elem.attr("id");
    var parent_id = parent_elem.attr("id");
    var child_id = child_elem.attr("id");

    $.ajax('/recipes/drop/',
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'source_id': source_id,
                'parent_id': parent_id,
                'child_id':c
            },
            success: function (data) {

                if (data.status === 'error') {
                    popup_message(parent_elem, data.msg, data.status, 2000);
                } else {
                    //alert(data.status);
                    source_elem.transition('zoom');
                    window.location.reload();
                    popup_message(parent_elem, "Moved Post", 'success', 2000);
                }
            },
            error: function (xhr, status, text) {
                error_message(parent_elem, xhr, status, text);
            }
        });

}

function drag_and_drop() {

    $(".droppable").droppable(
        {
            accept: ".recipes .item, .projects .item",
            drop: function (event, ui) {

                // Source post being dragged.
                var source = ui.draggable;

                // Parent post to drop into.
                var parent = $(this).closest(".recipes .item, .projects .item");
                if (!parent.length){
                    parent = $(this)
                }
                // Move target post to parent.
                move_object(parent, source)
            },
        });

    // Bind to any post object with the .draggable class
    $('.draggable').mousedown(function () {
        $(this).css('cursor', 'grabbing');
        var obj = $(this).closest('.recipes .item, .projects .item');

        obj.draggable(
        {
            addClasses: false,
            scroll: false,
            helper: 'clone',
            iframeFix: true,
            opacity: .7,
            containment: $('item'),
            revert:true,
            zIndex: 100,
            cursor:'grabbing'

        });
    });
}

