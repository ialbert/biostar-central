// Get cookies via jQuery.
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

// Get the csrftoken from the cookie.
csrftoken = getCookie('csrftoken');

function csrfSafeMethod(method) {
    // These HTTP methods do not require CSRF protection.
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}

$.ajaxSetup({
    // Ajax requests need to contain the csrf token.
    crossDomain: false, // obviates need for sameOrigin test
    beforeSend: function (xhr, settings) {
        if (!csrfSafeMethod(settings.type)) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});


// Create the comment panel.
function create_panel(elem, url) {
    // Remove the panel if it already exists.
    var panel = $('#panel')

    if (panel.length) {
        panel.remove();
    }
    var pk = elem.attr('data-pk')
    var container = elem.parent().parent()
    var panel = $('<div class="row" id="panel"></div>').load(url + pk + "/",
        function (response, status, xhr) {
            if (status == "error") {
                var msg = "Sorry but there was an error: ";
                alert(msg + xhr.status + " " + xhr.statusText);
            } else {
                $('#focus').focus();
            };
        }
    );
    container.after(panel);
}

function add_comment(elem) {
    return create_panel(elem, "/site/x/add_comment/")
}

function post_moderate(elem) {
    return create_panel(elem, "/site/x/post_moderate/")
}


function toggle_state(elem, vote_type) {
    // Toggles the state of the buttons and updates the label
    if (elem.hasClass('off')) {
        elem.removeClass('off');
        change = 1
    } else {
        elem.addClass('off');
        change = -1
    }
    // Counts stored in a sibling with class=count.
    var value = parseInt(elem.siblings('.count').text()) || 0
    elem.siblings('.count').text(value + change)
}

function pop_over(elem, msg, cls) {
    // A message panel that pops over the element.
    var div = '<div></div>'
    var tag = $(div).insertAfter(elem)
    tag.addClass('popover ' + cls)
    tag.text(msg)
    tag.delay(1000).fadeOut(1000, function () {
        $(this).remove()
    });
}

function vote_handler(elem, post_id, vote_type) {
    // Pre-emptitively toggle the button to provide feedback
    toggle_state(elem, vote_type)

    $.ajax('/site/x/vote/', {
        type: 'POST',
        dataType: 'json',
        data: {post_id: post_id, vote_type: vote_type},
        success: function (data) {
            if (data.status == 'error') { // Soft failure, like not logged in
                pop_over(elem, data.msg, data.status) // Display popover only if there was an error
                toggle_state(elem, vote_type) // Untoggle the button if there was an error
            } else {
                // Uncomment the line below to produce a popover on success as well.
                //pop_over(elem, data.msg, data.status)
            }

        },
        error: function () { // Hard failure, like network error
            pop_over(elem, 'Unable to submit vote!', 'error');
            toggle_state(elem, vote_type);
        }
    });
}

$(document).ready(function () {

    // Set focus on search field.
    var search = $('#search');

    // Move cursor at the end of the text in the field.
    if (search.length) {
        size = search.val().length;
        search.focus();
        search[0].setSelectionRange(size, size);
    }

    // Draw pie charts for users.
    $(".line").peity("line")

    // Initialize the sytnax highlighter.
    hljs.initHighlightingOnLoad();

    // Initialize the PageDown editor initializer
    var wmd = $('#wmd-input')
    if (wmd.length) {
        var converter = new Markdown.Converter();
        var editor = new Markdown.Editor(converter);
        editor.run();
    }

    // Register vote submission function.
    $('.vote').each(function () {
        $($(this)).click(function () {
            var elem = $(this);
            var post_id = elem.attr('data-post_id');
            var vote_type = elem.attr('data-type');
            vote_handler(elem, post_id, vote_type);
        });
    });

    // Authenticated user actions.
    $('.reply').each(function () {
        $(this).click(function (event) {
            event.preventDefault();
            add_comment($(this));
        });
    });

    $('.post_moderate').each(function () {
        $(this).click(function (event) {
            event.preventDefault();
            post_moderate($(this));
        });
    });

});
