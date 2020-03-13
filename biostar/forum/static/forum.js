function captcha() {

    var key = $("#recaptcha_key").val() || '';
    var captchaWidgetId = grecaptcha.render('captcha', {
        'sitekey': key,
        'theme': 'light'
    });

    return grecaptcha.getResponse(captchaWidgetId);
}


function init_pagedown(){
    var converter = new Markdown.getSanitizingConverter();
    var editor = new Markdown.Editor(converter);
    editor.run();
    return editor
}

function apply_vote(elem, post_uid, vote_type) {

    // Toggle the button to provide feedback
    elem.toggleClass("on");

    // The vote handler.
    vote_url = "/ajax/vote/";

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

function highlight_search(target, content_elem, stop_list) {

    // Find the target in the content.

    var target_list = target.replace(/\s{2,}/g, ' ').split(" ");

    // filter stoplist from target list.

    $.each(target_list, function (index, value) {

        var html = content_elem.html();
        var insert = "<span class='search-highlight'>" + value + "</span>";

        var new_html = html.replace(new RegExp(value, "ig"), insert);

        // Don't highlight any stop words.
        if ($.inArray(value, stop_list) === -1 && value.length >= 3) {
            console.log(value, insert, html);
            content_elem.html(new_html);
        }

    });
}

function highlight(text) {
    var con = markdownit({
        // ESCAPES when html=true
        html: false,
        highlight: function (str, lang) {
            if (lang && hljs.getLanguage(lang)) {
                try {
                    return '<pre class="hljs"><code>' +
                        hljs.highlight(lang, str, true).value +
                        '</code></pre>';
                } catch (__) {
                }
            }
            return '<pre class="hljs"><code>' + con.utils.escapeHtml(str) + '</code></pre>';
        }
    });

    var res = con.render(text);

    //preview.html(res);
    res.find('pre').addClass('language-bash');
    res.find('code').addClass('language-bash');

    //Prism.highlightAll()
}


function mark_spam(post_id, elem) {

    $.ajax('/ajax/report/spam/' + post_id + "/",
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: {},
            success: function (data) {
                //alert(elem.html());
                if (data.status === 'error') {
                    popup_message(elem.parent().parent(), data.msg, data.status);

                } else {
                    popup_message(elem.parent().parent(), data.msg, data.status);
                    $('#' + post_id).removeClass('open').addClass('spam');
                }

            },
            error: function (xhr, status, text) {
                error_message($('#' + post_id), xhr, status, text)
            }

        });


}


function moderate(elem, url) {

    event.preventDefault();
    //var elem = $(this);
    $('#modpanel').remove();

    // Could be a user or post uid
    var data_uid = elem.attr('data-value');

    var container = $('.moderate-insert[data-value="' + data_uid + '"]');
    var mod_url = url + data_uid + '/';
    //alert("FOOO");
    var page = $('<div id="modpanel"></div>').load(mod_url);
    container.after(page);
    //alert(container.html());
}

function similar_posts(elem) {
    var uid = elem.attr('post_uid');
    // Construct the similar posts link.
    var feed_url = '/similar/posts/' + uid + '/';
    var dimm_elem = $('#dim-similar');
    dimm_elem.dimmer('show');

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
                    dimm_elem.dimmer('hide');
                } else {
                    // Populate the feed.
                    dimm_elem.dimmer('hide');
                    elem.html(data.html);

                }
            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}

function change_subs(elem, value, $item) {
    var post_id = elem.closest('.post').data("value");
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


function activate_prism(elem){

    elem = elem || $(document);
    //if (!elem){
    elem.each('pre').addClass('language-bash');
    elem.each('code').addClass('language-bash');
    Prism.highlightAll();
}

function tags_dropdown(){

    $('.tag-field').dropdown({
        allowAdditions: true,
        // Get form field to add to
        onChange: function (value, text, $selectedItem) {
            // Get form field to add to
            var tagid = $(this).children('select').attr('field_id');
            var tag_field = $('#{0}'.f(tagid));
            // Add selected tag to field
            tag_field.val(value);

        }
    });

    $('.tag-field > input.search').keydown(function (event) {
        // Prevent submitting form when adding tag by pressing ENTER.
        if (event.keyCode === 13) {
            event.preventDefault();
        }
        // Set value with SPACE bar
        if (event.keyCode === 32) {
            event.preventDefault();
            //alert(alert( $(this).parent('.tag-field').children('select').attr('class')));
            $(this).closest('select').dropdown('set selected', $(this).val().trim());
            $(this).val('');
        }

    });
}

$(document).ready(function () {

    $('.mark-spam.item').click(function (event) {
        var post_id = $(this).closest('.post').attr('id');
        //alert($(this).closest('.post').attr('id'));
        mark_spam(post_id, $(this));
    });

    $('#similar-feed').each(function () {
        var elem = $(this);
        similar_posts(elem);
    });

    $('.ui.dropdown').dropdown();

    $(this).on('keyup', 'textarea', function (event) {
        var text = $(this).val();
        var highlighted = highlight(text);
        var form = $(this).closest('form');
        form.find('.preview').html(highlighted);
    });


    $(this).on('click', '#wmd-button-bar', function (event) {
        setTimeout(function () {
            var form = $(this).closest('form');
            var text = form.find('.textarea').val();
            var highlighted = highlight(text);
             //var form = $(this).closest('form');
             form.find('.preview').html(highlighted);
        }, 10);

    });

    $('#subscribe')
        .dropdown({
            action: 'hide',
            onChange: function (value, text, $item) {
                var elem = $(this);
                change_subs(elem, value, $item);
            }
        });



    $(".moderate-post").click(function (event) {
        event.preventDefault();
        var elem = $(this);
        moderate(elem, '/moderate/');
    });

    $(".moderate-user").click(function (event) {
        event.preventDefault();
        var elem = $(this);
        moderate(elem, '/accounts/moderate/');

    });

    $('.vote .upvote.button').popup({
            on: 'hover',
            content:'Upvote'
    });
    $('.vote .bookmark.button').popup({
            on: 'hover',
            content:'Bookmark '
    });
   $('.vote .accept.button').popup({
            on: 'hover',
            content:'Accept answer '
    });
    $('.voting .button').each(function (event) {
        var elem = $(this);
        var data_state = elem.attr('data-state');

        // Set the on class if the vote is selected.
        if (data_state === "1") {
            elem.addClass("on")
        }

        // Actions taken on vote click.
        $(this).click(function () {
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


    $('.ui.sticky').sticky();



    $('#show-answer').click(function () {
        $('.hidden-answer').toggle()
    });

    $('.trigger-highlight').each(function (event) {
        var elem = $(this);
        let container = $('#search-results');
        let query = container.data('query');
        let stop_words = container.data('stop');

        if (container.html() === undefined || container.html() === null) {
        } else {
            let stop_list = stop_words.split(',');
            highlight_search(query, elem, stop_list)
        }

    });

    tags_dropdown();
    activate_prism();
    init_pagedown();
    remove_trigger();
})
;
