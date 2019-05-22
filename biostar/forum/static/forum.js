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


function toggle_icon(elem) {

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


function mod_votecount(elem, k) {

    count = parseInt(elem.siblings('.count').text()) || 0;
    count += k;
    elem.siblings('.count').text(count)

};

function ajax_vote(elem, post_uid, vote_type, vote_url) {

    // Toggle the button to provide feedback
    toggle_icon(elem);

    // The vote handler.
    vote_url = "/vote/"

    $.ajax(vote_url, {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'post_uid': post_uid,
            'vote_type': vote_type,
        },

        success: function (data) {
            if (data.status === 'error') {
                toggle_icon(elem); // Untoggle the button if there was an error
                pop_over(elem, data.msg, data.status);
            } else {
                //pop_over(elem, data.msg, data.status)
            }

        },
        error: function () {
            toggle_icon(elem, vote_type);
        }
    });
}

function add_reply(elem) {

    // Remove body if it exists.
    $("#reply-row").remove();

    var msg_uid = elem.attr('data-value');
    var container = $("#reply-container-" + msg_uid);
    var reply_url = elem.attr("reply-url");

    var page = $('<div id="reply-row"></div>').load(reply_url);
    container.after(page);

}


function add_comment(elem) {

    var post_uid = elem.attr('data-value');
    var target = count_elem(post_uid)
    var url = "/create/comment/" + post_uid + "/"

    // remove comment body if exists.
    $("#comment-row").remove();


    var container = $("#comment-container-" + post_uid);
    var comment_url = elem.attr("comment-url");


    $.ajax(url, {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',
        success: function (data) {
            if (data.status === 'error') {
                // Untoggle the button if there was an error
                pop_message(target, data.msg, data.status);
                alert(data.msg);
            } else {
                pop_message(target, data.msg, data.status);
            }

        },
        error: function (xhr, status, text) {
            pop_message(target, text, "error")
        }
    });

    container.after(page);

};


function pop_message(elem, msg, cls) {
    var text = '<div></div>'
    var tag = $(text).insertBefore(elem)
    tag.addClass('popover ' + cls)
    tag.text(msg)
    tag.delay(1000).fadeOut(500, function () {
        $(this).remove()
    });
}


function count_elem(post_uid) {
    // The DOM element that stores the count
    elem = $("#count-" + post_uid)
    return elem;
}

function toggle_class(elem, post_uid) {

    var icon = elem

    // Toggles the state of the buttons and updates the label messages
    if (icon.hasClass('on')) {
        icon.removeClass("on");
        change = -1
    } else {
        icon.addClass("on")
        change = 1

    }
    // Increment the post score counter
    elem = count_elem(post_uid)
    var value = (parseInt(elem.text()) || 0) + change;
    elem.text(value)

};

function apply_vote(elem, post_uid, vote_type) {

    // Toggle the button to provide feedback
    toggle_class(elem, post_uid);

    // Message target.
    var target = count_elem(post_uid)

    // The vote handler.
    vote_url = "/ajax/vote/"

    $.ajax(vote_url, {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'post_uid': post_uid,
            'vote_type': vote_type,
        },

        success: function (data) {
            if (data.status === 'error') {
                // Untoggle the button if there was an error
                toggle_class(elem, post_uid);
                pop_message(target, data.msg, data.status);
            } else {
                pop_message(target, data.msg, data.status);
            }

        },
        error: function (xhr, status, text) {
            toggle_class(elem, post_uid);
            pop_message(target, text, "error")
        }
    });
}

$(document).ready(function () {

    $('select')
        .dropdown()
    ;

    $(".add-comment").click(function (event) {
        event.preventDefault();
        add_comment($(this));
    });

    $(".reply-button").click(function () {
        event.preventDefault();
        add_reply($(this));
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

    $("#handle").click(function () {
        var elem = $(this);
        //var uid = elem.attr("uid");
        var divHtml = elem.html(); //select's the contents of div immediately previous to the button

        // Textarea for user input
        var editableText = $("<div>" +
                             "<textarea name='handle' class='handle-text' rows='1'/>" +
                             "<span class='ui submit xsmall circular label'><i class='write icon'></i>Edit</span>" +
                             "</div>"
                            );
        //
        editableText.find('textarea[name="handle"]').val(divHtml);

        $.ajax()

        $(this).replaceWith(editableText);
        });

    $('.vote').each(function (event) {

        var elem = $(this);
        var vote_type = elem.attr('data-type');
        var data_state = elem.attr('data-state');

        // Popup information
        if (vote_type == "bookmark") {
            content = (data_state == "0") ? "Click to add bookmark" : "Click to remove bookmark"
        } else {
            content = (data_state == "0") ? "Click to add vote" : "Click to remove vote"
        }

        // Set the on class if the vote is selected.
        if (data_state == "1") {
            elem.addClass("on")
        }
        // Initialize each popup on the page.
        elem.popup({
                on: 'hover',
                content: content,
                delay: {
                    show: 800,
                    hide: 200
                }
            }
        );

        // Actions taken on vote click.
        $($(this)).click(function () {
            var elem = $(this);
            var post_uid = elem.attr('data-value');
            var data_type = elem.attr('data-type');
            apply_vote(elem, post_uid, data_type);

        });
    });


})
;
