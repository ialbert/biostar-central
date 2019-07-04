// As seen in: https://stackoverflow.com/questions/1038746/equivalent-of-string-format-in-jquery
// Usage: 'Hello {0} and {1}, and {0}'.f("World", "John");

String.prototype.format = String.prototype.f = function () {
    var s = this,
        i = arguments.length;

    while (i--) {
        s = s.replace(new RegExp('\\{' + i + '\\}', 'gm'), arguments[i]);
    }
    return s;
};

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


// Adds a comment to the post
function add_comment(elem) {

    var post_uid = elem.attr('data-value');
    var container = $("#comment-insert-" + post_uid);
    var url = "/new/comment/" + post_uid + "/"

    // Check for existing comment.
    var comment = $("#new-comment")

    if (comment.length) {
        // Remove comment if exists.
        comment.remove();
    } else {
        // Create a new comment.
        comment = $('<div id="new-comment"></div>')
    }

    // Insert into the page.
    container.after(comment);

    // Checks the size of the comment.
    function textarea_size_check() {
        var input = $("#comment-input");
        var size = input.val().length;
        if (size < 10) {
            popup_message(input, "More than 10 characters please!", "error");
        } else {
            //new_comment(post_uid, content);
            $("#comment-form").submit()
        }
    }

    // Submit form with CTRL-ENTER
    comment.keydown(function (e) {
        if ((e.ctrlKey || e.metaKey) && (e.keyCode == 13 || e.keyCode == 10)) {
            textarea_size_check()
        };
    });

    // Replace comment form from server
    comment.load(url, function (response, status, xhr) {
        if (status === 'success') {
            // Focus on the input
            $("#comment-input").focus();
            // Apply size check to submit button.
            $("#comment-submit").click(function (e) {
                e.preventDefault();
                textarea_size_check()
            });
        } else {
            error_message(elem, xhr, status, "")
            comment.remove()
        }
    });

};


function popup_message(elem, msg, cls, timeout) {
    timeout = typeof timeout !== 'undefined' ? timeout : 1000;
    var text = '<div></div>'
    var tag = $(text).insertBefore(elem)
    tag.addClass('popover ' + cls)
    tag.text(msg)
    tag.delay(timeout).fadeOut(500, function () {
        $(this).remove()
    });
}

// Triggered on network errors.
function error_message(elem, xhr, status, text) {
    popup_message(elem, "Error! readyState=" + xhr.readyState + " status=" + status + " text=" + text, "error", timeout = 5000)
}

function apply_vote(elem, post_uid, vote_type) {

    // Toggle the button to provide feedback
    elem.toggleClass("on")

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
                elem.toggleClass("on")
                popup_message(elem, data.msg, data.status);
            } else {
                // Success
                //popup_message(elem, data.msg, data.status);
                // Increment the post score counter
                var score = $("#score-" + post_uid)
                var value = (parseInt(score.text()) || 0) + parseInt(data.change) || 0;
                score.text(value)
            }

        },
        error: function (xhr, status, text) {
            elem.toggleClass("on")
            error_message(elem, xhr, status, text)
        }
    });
}

function remove_trigger() {
    // Makes site messages dissapear.
    $('.remove').delay(2000).slideUp(800, function () {
        $(this).remove();
    });
}

function edit_post(post_uid) {
    var edit_url = '/ajax/edit/';
    var form_elem = $('#inplace-form-' + post_uid);
    var edited = form_elem.find('textarea').val();
    $.ajax(edit_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'post_uid': post_uid,
                'content': edited
            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(form_elem, data.msg, data.status, 3000);
                } else {
                    // Hide form
                    $('#inplace-form-' + post_uid).hide();
                    // Replace with edited data
                    $('#inplace-' + post_uid).html(data.msg).show().focus();
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

$(document).ready(function () {

    $('.ui.dropdown').dropdown();

    $('.tag-field').dropdown({
        allowAdditions: true,
        onChange: function(value, text, $selectedItem) {
            // Get form field to add to
            var tagid = $("#tag-menu").attr('field_id');
            var tag_field = $('#{0}'.f(tagid));
            // Add selected tag to field
            //alert(value);
            tag_field.val(value);
    }
    });
    $('.tag-field >input.search').keydown(function(event) {
        // Prevent submitting form when adding tag by pressing ENTER.
        if (event.keyCode === 13){
            event.preventDefault();
        }
        // Set value with SPACE bar
        if (event.keyCode === 32){
            event.preventDefault();
            $("#tag-menu").dropdown('set selected', $(this).val().trim());
            $(this).val('')
        }

    });

    $('.inplace').click(function () {
        // Hide content
        $(this).hide();
        var post_uid = $(this).attr('post_uid');
        // Exposes form
        $('#inplace-form-' + post_uid).show().focus();
    });

    $('.inplace-form button.cancel').click(function (event) {
        event.preventDefault();
        var post_uid = $(this).closest('.inplace-form').attr('post_uid');
        // Hide form and show content
        $('#inplace-form-' + post_uid).hide();
        $('#inplace-' + post_uid).show().focus();
    });

    $('.inplace-form textarea').keyup(function (event) {
        var post_uid = $(this).closest('.inplace-form').attr('post_uid');

        // Submit edit when pressing CTRL-ENTER
        if ((event.ctrlKey || event.metaKey) && (event.keyCode === 13 || event.keyCode === 10)) {
            event.preventDefault();
            edit_post(post_uid);
            return;
        }
        // Leave when pressing ESC
        if (event.keyCode === 27){
            // Hide form and show content
            $('#inplace-form-' + post_uid).hide();
            $('#inplace-' + post_uid).show().focus();

        }

    });

    $('.inplace-form button.save').click(function () {
        // Submit edit when clicking save
        event.preventDefault();
        var post_uid = $(this).closest('.inplace-form').attr('post_uid');
        edit_post(post_uid);
    });

    $('.inplace-edit').click(function () {
        var post_uid = $(this).attr('post_uid');
        $('#inplace-' + post_uid).hide();
        $('#inplace-form-' + post_uid).show().focus();
    });

    $('#subscribe')
        .dropdown({
          action:'hide',
          onChange: function (value, text, $item) {
          var elem = $(this);

          console.log(elem.id);

          // Get the root id
          var post_id = elem.attr("data-uid");
          // Currently selected item
          var active = $('#active');
          // Subscription url
          var subs_url = '/ajax/subscribe/';
          $.ajax(subs_url,
              {
                  type: 'POST',
                  dataType: 'json',
                  ContentType: 'application/json',
                  data: {
                      'root_uid': post_id,
                      'sub_type': value
                  },
                  success: function (data) {
                      if (data.status === 'error') {
                          popup_message(elem, data.msg, data.status);
                      } else {
                          // Replace current item with the select one.
                          active.text($item.text());
                      }

                      },
                  error: function (xhr, status, text) {
                    error_message(elem, xhr, status, text)
                  }
              })
          }
        });

    $(".add-comment").click(function (event) {
        event.preventDefault();
        add_comment($(this));
    });


    remove_trigger();

    $('#search').keyup(function (event) {
        var search_url='/ajax/search/';
        var query = $(this).val();

        $.ajax(search_url, {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'query':query
            },

            success: function (data) {
                if (data.status === 'error') {
                    // Untoggle the button if there was an error
                } else {
                    // Success
                    //popup_message(elem, data.msg, data.status);
                    // Increment the post score counter
                    var res = $('#results');
                    //alert(data.html);
                    res.html(data.html)
                }

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        });


    });

    $(".moderate-post").click(function (event) {
        event.preventDefault();
        var elem = $(this);

        $('#modpanel').remove();

        // Could be a user or post uid
        var data_uid = elem.attr('data-value');

        var container = $("#moderate-insert-" + data_uid);
        var mod_url = '/moderate/' + data_uid + '/';

        var page = $('<div id="modpanel"></div>').load(mod_url);
        container.after(page)

    });

    $(".moderate-user").click(function (event) {
        event.preventDefault();
        var elem = $(this);

        $('#modpanel').remove();

        // Could be a user or post uid
        var data_uid = elem.attr('data-value');

        var container = $("#moderate-insert-" + data_uid);
        var mod_url = '/accounts/moderate/'+ data_uid + '/';

        var page = $('<div id="modpanel"></div>').load(mod_url);
        container.after(page)
    });

    $('.vote').each(function (event) {
        var elem = $(this);
        var data_state = elem.attr('data-state');

        // Set the on class if the vote is selected.
        if (data_state === "1") {
            elem.addClass("on")
        }

        // Actions taken on vote click.
        $($(this)).click(function () {
            var elem = $(this);
            var post_uid = elem.attr('data-value');
            var data_type = elem.attr('data-type');
            apply_vote(elem, post_uid, data_type);
        });
    });

    $("#form-errors .error").each(function () {

        var elem = $(this);
        // Get errored out field id and label
        var field_id = elem.attr('data-value');
        var field_label = elem.attr('label');
        // Get the error message
        var message = elem.attr("message");
        // Select field in the form using it's id
        var field = $(field_id);
        // Add an 'error' to '.ui.field' to turn it red.
        field.closest(".field").addClass("error");
        // Insert the error message
        field.before('<div class="ui small red message"> {1}</div>'.f(field_label, message))
    });


    $('pre').addClass('language-bash');
        Prism.highlightAll();

})
;
