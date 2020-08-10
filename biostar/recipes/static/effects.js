


function move_object(parent_elem, source_elem, next_elem, url){

    var source_id = source_elem.attr("id");
    var parent_id = parent_elem.attr("id");
    var next_id = next_elem.attr("id");

    $.ajax(url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'source_id': source_id,
                'parent_id': parent_id,
                'next_id': next_id,
            },
            success: function (data) {

                if (data.status === 'error') {
                    popup_message(parent_elem, data.msg, data.status, 2000);
                } else {
                    //alert(data.status);
                    //window.location.reload();
                    parent_elem.after(source_elem);
                    //popup_message(parent_elem, "Moved Object", 'success', 2000);
                    source_elem.transition('pulse').transition('glow')

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
            accept: ".recipe.item, .projects .item",
            drop: function (event, ui) {

                // Source object being dragged.
                var source = ui.draggable;

                // Get the object to drop under
                var parent = $(this).closest(".recipe.item, .projects .item");
                if (!parent.length){
                    parent = $(this)
                }

                // Get the object that will be below.
                var next = parent.next();

                // Resolve the url.
                if (parent.closest(".projects").hasClass('projects') || source.closest(".projects").hasClass('projects')){
                    url = '/project/drop/'
                }else{
                    url = '/recipe/drop/';

                }

                if (!next.attr("id")){
                    next = parent.next().next();
                }

                move_object(parent, source, next, url)
            },
        });

    // Bind to any object with the .draggable class
    $('.draggable').mousedown(function () {
        $(this).css('cursor', 'grabbing');
        var obj = $(this).closest('.recipe.item, .projects .item');

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

