// using jQuery
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie != '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) == (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}

csrftoken = getCookie('csrftoken');

function csrfSafeMethod(method) {
    // these HTTP methods do not require CSRF protection
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}

$.ajaxSetup({
    crossDomain: false, // obviates need for sameOrigin test
    beforeSend: function(xhr, settings) {
        if (!csrfSafeMethod(settings.type)) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});

// comment add
function user_comment_click(elem) {

    // remove comment body if exists.
    $("#comment-form").remove();

    var post_id = elem.attr('data-value')
    var container = elem.closest("table")

    var csrf_html = $('#csrf_token').find("input[name='csrfmiddlewaretoken']").parent().html()

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

function toggle_button(elem) {
    // Toggles the state of the buttons and updates the label messages
    if (elem.hasClass('off')) {
        elem.removeClass('off');
    } else {
        elem.addClass('off');
    }
}

function pop_over(elem, msg, cls) {

    elem.append('<div></div>')
    var tag = elem.children('div').last()
    tag.addClass('vote-popover ' + cls)
    tag.text(msg)
    tag.delay(1000).fadeOut(1000, function () {
        $(this).remove()
    });
}

function ajax_vote(elem, post, type) {
    // Pre-emptitively toggle the button to provide feedback
    toggle_button(elem)

    $.ajax('/x/vote/', {
        type: 'POST',
        dataType: 'json',
        data: {post: post, type: type},
        success: function (data) {
            if (data.status == 'error') { // Soft failure, like not logged in
                pop_over(elem, data.msg, data.status) // Display popover only if there was an error
                toggle_button(elem) // Untoggle the button if there was an error
            } else {
                //pop_over(elem, "Vote Success", "success")
            }

        },
        error: function () { // Hard failure, like network error
            pop_over(elem, 'Unable to submit vote!', 'error');
            toggle_button(elem);
        }
    });
}

$(document).ready(function () {
    var tooltip_options = {};

    // This detects the user id
    var user_id = $("#user_id").val()

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

    $('.vote').each(function () {

        $($(this)).click(function () {
            var elem = $(this)
            var post = 100
            var type = elem.attr('data-type')
            ajax_vote(elem, post, type);
        });
    });


});
