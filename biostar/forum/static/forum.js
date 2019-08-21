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


function edit_post(uid) {

    var edit_url = '/ajax/edit/' + uid + '/';
    // Rendered form element
    var form_elem = $('.edit-form[data-value="'+ uid +'"]');
    // Inplace form container
    var form_container = $('inplace[data-value="'+ uid +'"]');
    // Hidden elements
    var hidden =  $('.hide-on-edit[data-value="'+ uid +'"]');

    // Post title inside of the form
    var title = $('.title[data-value="'+ uid +'"]');
    var content = $('.content[data-value="'+ uid +'"]');
    var post_type = $('#inplace-type').dropdown('get value');
    var tag_val = $('.tag-field').dropdown('get value');

    // Current post content and title to replace
    // with returned values.
    var post_content = $('.editable[data-value="'+ uid +'"]');
    var post_title = $('.post-title[data-value="'+ uid +'"]');
    var post_tags = $('.post-tags[data-value="'+ uid +'"]');

    //alert(content.val());
    // Title is null and type are null
    // meaning current post is not top level.
    if (title.val() == null){
        title.val('')
    }
    if (!($.isNumeric(post_type))){
        post_type = -1
    }

    $.ajax(edit_url,
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            traditional: true,
            data: {
                'content': content.val(),
                'title':title.val(),
                'type':post_type,
                'tag_val':tag_val,

            },
            success: function (data) {
                if (data.status === 'error') {
                    popup_message(form_elem, data.msg, data.status, 3000);
                } else {
                    // Clear and hide inplace form
                    form_elem.html('');
                    form_container.hide();

                    // Show hidden items
                    hidden.show();

                    // Replace current post info with edited data
                    post_content.html(data.html).show().focus();
                    post_title.html(data.title).show();
                    post_tags.html(data.tag_html).show();

                    // Highlight the markdown in content
                    post_content.find('pre').addClass('language-bash');
                    post_content.find('code').addClass('language-bash');
                    Prism.highlightAll();
                }
            },
            error: function (xhr, status, text) {
                error_message(form_elem, xhr, status, text)
            }
        })
}

function cancel_inplace(uid){

    var inplace_content = $('inplace[data-value="'+ uid +'"]');
    //var inplace_title = $('inplace-title[data-value="'+ uid +'"]');

    //var title = $('.editable-title[data-value="'+ uid +'"]');
    var content = $('.editable[data-value="'+ uid +'"]');
    //var content = $('#content-' + uid);
    var hidden = $('.hide-on-edit[data-value="'+ uid +'"]');
    //Delete the form
    inplace_content.html("");
    //inplace_title.html("");
    // Hide the container
    // Show original content
    content.show();
    //Show any blocked element
    hidden.show();

}


function inplace_post(elem){

    var uid = elem.data("value");
    var hidden =  $('.hide-on-edit[data-value="'+ uid +'"]');
    var form_container = $('inplace[data-value="'+ uid +'"]');
    var url = '/inplace/post/' + uid +'/';

    $.ajax(url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            success: function (data) {
                if (data.status === 'error') {
                    alert(data.status);
                    alert(data.msg);
                    popup_message(elem, data.msg, data.status, 3000);
                } else {
                    elem.hide();
                    hidden.hide();
                    form_container.html(data.inplace_form);
                    form_container.show().find('textarea').focus();
                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}


function  search(query, elem, search_url) {
    var res = $('#results');
    var container = $('#contain-search');

    container.addClass('loading search');
    res.width(container.width());
    res.addClass('ui search message');
    res.html('Searching ...');

    $.ajax(search_url, {
        type: 'GET',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'query': query
        },

        success: function (data) {
            if (data.status === 'error') {
                popup_message($(this), data.msg, data.status);
            } else {
                // Success
                //alert(data.html);
                res.removeClass('ui message');
                res.html(data.html);
                container.removeClass('loading search');
            }

        },
        error: function (xhr, status, text) {
            error_message(elem, xhr, status, text)
        }
    });
}


$(document).ready(function () {

    $(this).click(function(event) {
        var res = $('#results');
        if(typeof res.html() === 'undefined'){
            return
        }
        if (res.html().length > 0){
            res.html('');
            res.removeClass('ui message');
        }
    });

    $('#similar-feed').each(function () {
        var elem = $(this);
        //elem.html('asshole');
        // Fetch feed from url url
        var uid = elem.attr('post_uid');
        var feed_url = '/similar/posts/' + uid + '/';


        $.ajax(feed_url,
            {
                type: 'GET',
                dataType: 'json',
                ContentType: 'application/json',
                data: {
                    'uid': uid
                },
                success: function (data) {
                    if (data.status === 'error') {
                        popup_message(elem, data.msg, data.status);
                    } else {
                        // Populate the feed.
                        elem.html(data.html)
                    }
                },
                error: function (xhr, status, text) {
                    error_message(elem, xhr, status, text)
                }
            })
    });

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

    $('.editable').click(function (event) {
         if (event.metaKey || event.ctrlKey){
             inplace_post($(this))
         }
    }).dblclick(function (event) {
         inplace_post($(this))
    });

    $('.inplace-click').click(function (event) {
        var uid = $(this).data('value');
        var elem = $('.editable[data-value="'+ uid +'"]');
        inplace_post(elem);
    });

    $(this).on('keyup', '.edit-form textarea', function (event) {

        var uid = $(this).data('value');
        // Submit form with CTRL-ENTER
        if (event.ctrlKey && (event.keyCode === 13 || event.keyCode === 10)) {
            edit_post(uid);
            return
        }

        var md = markdownit();
        var text = $(this).val();
        var html_preview = md.render(text);
        //var html_preview = Prism.highlight(md.render(text), Prism.languages.bash, 'language-bash');
        var html_container = $('#html-preview-'+ uid);

        html_container.html(html_preview);
        //alert("test");
        html_container.find('pre').addClass('language-bash');
        html_container.find('code').addClass('language-bash');
        Prism.highlightAll();



    });

    $(this).keyup(function (event) {
        if (event.keyCode === 27){
            $('.inplace').each(function () {
                event.preventDefault();
                var uid = $(this).data("value");
                cancel_inplace(uid);
            });
        }

    });

    $(this).on('click', '.edit-form .cancel', function(){
        event.preventDefault();
        var uid = $(this).data("value");
        cancel_inplace(uid);
     });

    $(this).on('click', '.edit-form .save', function(){
        var uid = $(this).data("value");
        event.preventDefault();
        edit_post(uid);
     });


    $('#digest').dropdown({
          action:'hide',
          onChange: function (value, text, $item) {
          var elem = $(this);

          console.log(elem.id);

          // Get the root id
          var uid = elem.data("uid");
          // Currently selected item
          var active = $('#digest-active');
          var icon_container = $('#digest-icon');
          var icon_str = $item.data('icon');
          // Subscription url
          var digest_url = '/ajax/digest/' + uid + '/';
          $.ajax(digest_url,
              {
                  type: 'POST',
                  dataType: 'json',
                  ContentType: 'application/json',
                  data: {
                      'pref': value
                  },
                  success: function (data) {
                      if (data.status === 'error') {
                          popup_message(elem, data.msg, data.status);
                      } else {
                          // Replace current item with the select one.
                          active.text($item.text());
                          icon_container.removeClass();
                          icon_container.addClass(icon_str);
                      }
                      },
                  error: function (xhr, status, text) {
                    error_message(elem, xhr, status, text)
                  }
              })
          }
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
          var active = $('#sub-active');
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
        var query = $(this).val();
        var search_url = $(this).attr('url');
        // Only preform searches when pressing ENTER
        if (event.keyCode === 13){
            search(query, $(this), search_url);
        }
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

    $('a').click(function () {
        $(this).transition('pulse', 295);
    });

    $(this).on('click', '.show-preview', function() {
        var uid = $(this).data('value');
        var preview = $('.preview-'+uid);
        preview.transition('slide down', 400);
        preview.find('pre').addClass('language-bash');
        preview.find('code').addClass('language-bash');
        Prism.highlightAll();

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

    $('.display-answer').click(function() {
        $('.answer-form').transition('slide down');
    });

    $('pre').addClass('language-bash');
    $('code').addClass('language-bash');
    Prism.highlightAll()
})
;
