function captcha() {

    var key = $("#recaptcha_key").val();

    if (key) {
        var captchaWidgetId = grecaptcha.render('captcha', {
            'sitekey': key,
            'theme': 'light'
        });
        return grecaptcha.getResponse(captchaWidgetId);

    }
}


function apply_vote(vote_elem) {

    var post = vote_elem.closest('.post');
    var uid = post.data("value");
    var type = vote_elem.data("value");
    var icon = vote_elem.find(".icon");

    $.ajax("/ajax/vote/", {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',
        data: {
            'post_uid': uid,
            'vote_type': type,
        },

        success: function (data) {
            if (data.status === 'error') {
                // Untoggle the button if there was an error
                icon.removeClass("on");
                popup_message(vote_elem, data.msg, data.status);
                return
            }
            // Success
            var score = vote_elem.closest('.voting').find(".score");
            var value = (parseInt(score.text()) || 0) + parseInt(data.change) || 0;
            score.text(value);
            icon.toggleClass("on");

        },
        error: function (xhr, status, text) {
            //icon.toggleClass("on");
            error_message(vote_elem, xhr, status, text)
        }
    });
}

function remove_trigger() {
    // Makes site messages dissapear.
    $('.remove').delay(2000).slideUp(800, function () {
        $(this).remove();
    });
}

function highlight(text) {

    var con = markdownit({
        // escape html in previews with html=false.
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

    return con.render(text);
    //preview.html(res);
}


function moderate(uid, container, url) {

    //event.preventDefault();
    //var elem = $(this);
    $('#modpanel').remove();

    // Could be a user or post uid
    var page = $('<div id="modpanel"></div>');
    container.after(page);
    page.hide();
    page.load(url, function (response, status, xhr) {
        page.show(300);
    });
    //page.show(1000)
    //alert("GGG")


}


function disable_emails(user_id, elem){

    var url = '/email/disable/{0}/'.format(user_id);
    $.ajax(url, {
        type: 'POST',
        dataType: 'json',
        ContentType: 'application/json',

        success: function (data) {
            if (data.status === 'error') {
                popup_message(elem, data.msg, data.status);
                return
            }
            // Success
            popup_message(elem, data.msg, data.status);
        },
        error: function (xhr, status, text) {
            //icon.toggleClass("on");
            error_message(elem, xhr, status, text)
        }
    });

}

function similar_posts(elem) {
    var uid = elem.attr('post_uid');
    // Construct the similar posts link.
    var feed_url = '/similar/posts/' + uid + '/';
    var dimm_elem = $('#similar');
    dimm_elem.dimmer('show');

    $.ajax(feed_url,
        {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
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

function change_subs(elem, value) {
    var post_id = elem.closest('.post').data("value");
    // Currently selected item
    // var active = $('#sub-active');
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
                    // active.text($item.text());
                }

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        })
}


function activate_prism(elem) {

    elem = elem || $(document);
    //if (!elem){
    elem.each('pre').addClass('language-bash');
    elem.each('code').addClass('language-bash');
    Prism.highlightAll();
}

function tags_dropdown() {

    $('.tags').dropdown({
        allowAdditions: true,
        // Get form field to add to
        onChange: function (value, text, $selectedItem) {
            // Get form field to add to
            var field = $(this).find("select").data("value");
            var tag_field = $('#{0}'.f(field));
            // Add selected tag to field
            // Set text instead of value
            value = $('<div/>').text(value).html();
            tag_field.val(value);
        }
    });
    $('.tags > input.search').keydown(function (event) {

        // Prevent submitting form when adding tag by pressing ENTER.
        var ek = event.keyCode || event.which;
        var value = $(this).val().trim();

        // Get a list of delimiters
        var delimiters = $('#field-tags').data('delimiters').split(',');

        if (delimiters.indexOf(String(ek)) !== -1) {
            // Escape the text before settings value.
            value = $('<div/>').text(value).html();
            event.preventDefault();
            $(this).closest('.tags').dropdown('set selected', value);
            $(this).val('');
            return value
        }
    })

}

function highligh_preview(form, text) {
    var highlighted = highlight(text);

    form.find('.preview').html(highlighted);
    form.find('pre').addClass('language-bash');
    form.find('code').addClass('language-bash');

    Prism.highlightAll();


}


$(document).ready(function () {

    $('#similar-feed').each(function () {
        var elem = $(this);
        similar_posts(elem);
    });

    $('.ui.dropdown').dropdown();

    $('form .preview').each(function () {
        var text = $(this).closest('form').find('.wmd-input').val();
        var form = $(this).closest('form');
        setTimeout(function () {
            highligh_preview(form, text);
        }, 150)
    });

    $(this).on('keyup', 'textarea', function (event) {
        var text = $(this).val();
        var form = $(this).closest('form');
        setTimeout(function () {
            highligh_preview(form, text);
        }, 50)
    });

    $(this).on('click', '#wmd-button-bar', function (event) {
        var form = $(this).closest('form');
        var text = form.find('textarea').val();
        highligh_preview(form, text);

    });

    $("#wmd-button-bar").click(function (event) {
        var form = $(this).closest('form');
        var text = form.find('textarea').val();
        highligh_preview(form, text);

    });

    $("body").on("click", '.pagedown-image-upload.show .submit-input', function () {

        var form = $(this).closest('form');
        setTimeout(
            function () {
                var text = form.find('textarea').val();
                highligh_preview(form, text);
            }, 500);
    });

    $('#subscribe').dropdown();

    $(this).unbind("click").on('click', '#subscribe .item', function (event) {
        var elem = $(this).closest('#subscribe');
        var value = $(this).data('value');
        change_subs(elem, value);

    });


    $(".profile .moderate").click(function (event) {
        event.preventDefault();
        var profile = $(this).closest('.profile');
        var uid = profile.data("value");
        var container = profile.find("#mod");
        var url = '/accounts/moderate/{0}/'.format(uid);
        moderate(uid, container, url)

    });
    $(".profile .disable-emails").click(function (event) {
        event.preventDefault();
        var profile = $(this).closest('.profile');
        var uid = profile.data("value");
        disable_emails(uid, profile)

    });
    $(".post .moderate").click(function (event) {
        event.preventDefault();
        var post = $(this).closest('.post');
        var uid = post.data("value");
        var container = $(this).closest('.post > .body > .content >.inplace');
        var url = '/moderate/{0}/'.format(uid);

        moderate(uid, container, url);
        cancel_inplace();


    });

    $("[data-value='upvote']").popup({
        on: 'hover',
        content: 'Upvote'
    });
    $(".draggable").popup({
        on: 'hover',
        content: 'Drag and Drop'
    });

    $("[data-value='bookmark']").popup({
        on: 'hover',
        content: 'Bookmark '
    });

    $("[data-value='accept']").popup({
        on: 'hover',
        content: 'Accept answer '
    });
    $('.voting button').each(function (event) {

        var elem = $(this);
        var data_state = elem.data('state');
        data_state = '{0}'.format(data_state);
        // Set the on class if the vote is selected.
        if (data_state === "1") {
            elem.find('.icon').addClass("on")
        }

    });
    $('.voting .button').unbind("click").click(function () {
        apply_vote($(this));
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

    tags_dropdown();

    $('pre').addClass('language-bash');
    $('code').addClass('language-bash');
    Prism.highlightAll();


})
;
