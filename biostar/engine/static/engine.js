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
    beforeSend: function (xhr, settings) {
        if (!csrfSafeMethod(settings.type)) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});

function mod_votecount(elem, k) {
    count = parseInt(elem.siblings('.count').text()) || 0
    count += k
    elem.siblings('.count').text(count)
}

function pop_over(elem, msg, cls) {
    var text = '<div></div>'
    var tag = $(text).insertAfter(elem)
    tag.addClass('ui message ' + cls)
    tag.text(msg)
    tag.delay(1000).fadeOut(1000, function () {
        $(this).remove()
    });
}

function submit_comment (elem) {

    var comment_form = $("#comment-form");
    var action_url = comment_form.attr("action");
    var container = $("#comment-row");

    $.ajax({
        url : action_url,
        type : "POST",
        ContentType: 'application/json',
        data : comment_form.serialize(),

        success: function (data) {
            if (data.status == 'error') {
                alert(data.msg)
            } else {
                //pop_over(elem, data.msg, data.status)
                location.reload(container);
            }

        },
        error: function () { // Hard failure, like network error
             alert(data.msg)
        }

    });
        return false;
};

// Triggered on moderation.
function moderate(elem) {

    $('#modpanel').remove();

    // Could be a user or post uid
    var data_uid = elem.attr('data-value');

    var container = $("#mod-container-" + data_uid);
    var mod_url = elem.attr('mod-url');

    var page = $('<div id="modpanel"></div>').load(mod_url);
    container.after(page)

};


function add_answer(elem) {

};

function toggle_button(elem) {

    var icon = elem.children("i.icon");

    // Toggles the state of the buttons and updates the label messages
    if (icon.hasClass('outline')) {
        icon.removeClass('outline');
        change = 1
    } else {
        icon.addClass('outline');
        change = -1
    }
    mod_votecount(elem, change)
}

function ajax_vote(elem, post_uid, vote_type, vote_url) {

    // Pre-emptitively toggle the button to provide feedback
    toggle_button(elem);

    $.ajax(vote_url, {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'post_uid': post_uid,
            'vote_type': vote_type,
            },

        success: function (data) {
            if (data.status == 'error') {
                toggle_button(elem); // Untoggle the button if there was an error
                pop_over(elem, data.msg, data.status);
            } else {
                //pop_over(elem, data.msg, data.status)
            }

        },
        error: function () {
            toggle_button(elem, vote_type);
        }
    });
}

function add_comment(elem) {

    // remove comment body if exists.
    $("#comment-row").remove();

    var post_uid = elem.attr('data-value');
    var project_uid = elem.attr('project-uid');
    var container = $("#comment-container-"+ post_uid);
    var comment_url = elem.attr("comment-url");
    var csrf_html = jQuery("[name=csrfmiddlewaretoken]").val();

    // Going to be refactored out and loaded separately
    container.after(`<div id="comment-row" class="ui basic segment inputcolor">
    <form id="comment-form" class="ui form" action=${comment_url}  method="post">
        <input type="hidden" name="parent_uid" id="parent_uid" value=${post_uid} />
        <input type="hidden" name="project_uid" id="project_uid" value=${project_uid} />
        <input type="hidden" name="csrfmiddlewaretoken" value=${csrf_html} />

        <div class="">
            <div id="wmd-button-bar-2"></div>
            <textarea class="wmd-input-2" id="wmd-input-2"  name="content" rows="6"></textarea>
        </div>
        <div>

            <a class="ui submit green button" onclick="return submit_comment($(this));">
                <i class="check icon"></i>Add Comment
            </a>
            <a class="ui orange right floated button" onclick="javascript:obj=$(\'#comment-row\').remove();">
            <i class="undo icon"></i> Cancel
            </a>
        </div>
    </form>
    </div>`
    );

    var converter = new Markdown.Converter();
    var editor = new Markdown.Editor(converter, '-2');
    editor.run();

};


$(document).ready(function () {


     $('select')
        .dropdown()
    ;

//    $(".items > .item").click(function (event) {
//        var obj = $(this).find("a:first");
//        if (typeof obj !== 'undefined') {
//            window.location = obj.attr("href");
//       }
//    });

    $(".copy-data").click(function (event) {
        event.preventDefault();
        var elem = $(this);
        var data_uid = elem.attr('data-uid');
        var copy_url = elem.attr('copy-url');

        $.ajax(copy_url, {
                type: 'GET',
                dataType: 'json',
                ContentType: 'application/json',
                data: {data_uid: data_uid},
                success: function (data) {
                $("#copy-message-"+ data_uid).append(`
                <div class="ui basic segment">
                    <div class="ui fluid green message">
                    ${data.message}
                    </div>
                </div>

                `).fadeOut(2000);
                        },
                });

    });

    $(".add-comment").click(function (event) {
        add_comment($(this));
        });

    $(".add-answer").click(function (event) {
        add_answer($(this));
        });

    $(".moderate-post").click(function (event) {
        moderate($(this));
        });

    $(".moderate-user").click(function (event) {
        moderate($(this));
        });


    $('.vote').each(function () {

        $($(this)).click(function () {
            var elem = $(this);
            var post_uid = elem.attr('data-post_uid');
            var vote_type = elem.attr('data-type');
            var vote_url = elem.attr('vote-url');

            ajax_vote(elem, post_uid, vote_type, vote_url);
        });
    });

});
