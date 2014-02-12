// comment add
function user_comment_click(elem) {

    // remove comment body if exists.
    $("#comment-form").remove();

    post_id = elem.attr('data-value')
    container = elem.closest("table")

    csrf_html = $('#csrf_token').find("input[name='csrfmiddlewaretoken']").parent().html()

    container.append('<tr><td colspan="2">\
    <form role="form" action="/p/new/comment/' + post_id + '/" method="post" id="comment-form" id="comment-form">' + csrf_html + '\
        <div class="form-group">\
        <textarea class="input-xlarge span8" id="comment-box" name="content" rows="3"></textarea></div> \
        <div><a class="btn btn-success" href=\'javascript:document.forms["comment-form"].submit()\'><i class="icon-comment"></i> Add comment</a>          \
        <a class="btn btn-warning pull-right" onclick="javascript:obj=$(\'#comment-form\').remove();"><i class="icon-remove"></i> Cancel</a>   </div>       \
    </form>            \
    </td></tr>'
    )
    CKEDITOR.replace('comment-box');
}

function anon_comment_click(elem) {
    container = elem.closest("table")
    elem.css("background-color", "red");
    $("#comment-box").remove();
    container.append('<tr id="comment-box">\
    <td colspan="2">\
        <div class="alert alert-warning">Please log in to comment</div>\
    </td></tr>'
    )

}

function anon_actions() {
    $('.add-comment').each(function () {
        elem = $(this)
        //console.log(elem)
        //callback function defined in /static/js/widgets.js
        elem.click(function () {
            show_add_comment($(this));
        });
    });
}


$(document).ready(function () {
    var tooltip_options = {};

    user_id = $("#user_id").val()


    // Register tooltips.
    $('.tip').tooltip(tooltip_options)

    // Register comment adding.
    if (user_id) {
        $('.add-comment').each(function () {
            $(this).click(function () {
                user_comment_click($(this));
            });
        });
    } else {
        $('.add-comment').each(function () {
            $(this).click(function () {
                anon_comment_click($(this));
            });
        });
    }

});
