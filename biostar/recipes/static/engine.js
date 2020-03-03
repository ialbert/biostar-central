//function


function update_job(data){
    if (data.state_changed) {
                    // TODO: Taking this type of styling out
                    $('.job-container-' + job_uid).html(data.html);
                    var job_item = $('.job-item-' + job_uid);
                    job_item.transition('pulse');

                }
                if (data.is_running) {
                    $('.loader[data-uid="' + job_uid + '"]').html('<div class="ui log message">\n' +
                        '                     <span class="ui active small inline loader"></span>\n' +
                        '                     <span>Running</span>\n' +
                        '</div>');
                    $('#job-link[data-uid="' + job_uid + '"]').attr("href", '/job/view/' + job_uid + '/#log')

                } else {
                    $('.loader[data-uid="' + job_uid + '"]').html("");
                    $('#job-link[data-uid="' + job_uid + '"]').attr("href", '/job/view/' + job_uid)
                }


                //alert(data.redir)
                imag.replaceWith(data.img_tmpl);

                if (data.stdout) {
                    var stdout = $('#stdout');
                    //alert($('#stdout').html());
                    stdout.text(data.stdout);
                    //alert(data.stdout);
                    //alert(stdout.html())
                }

                if (data.stderr) {
                    var stderr = $('#stderr');
                    stderr.text(data.stderr);

                }
                if (data.redir && $("#view").length) {
                    window.location.replace(data.redir + "#flist");
                    window.location.reload()
                }
}

function check_jobs() {

    // Look at each job with a 'check_back' tag
    $('.check_back').each(function () {

        var job = $(this).closest('.job');
        // Get the uid and state
        var uid = job.data('value');

        // Get the current job state.
        var state = $(this).data("state");

        // Get the element containing job image.
        var imag = job.find('#img:first');

        // Construct the url
        var url = '/ajax/check/job/{0}/'.format(uid);

        var link = job.find(".link");

        var stdout = $('#stdout pre:first');

        var loader = $('#stdout .loader');

        var stderr = $('#stderr');

        // Bail out when a job uid is not provided.
        if (uid === null || uid === undefined) {
            return
        }

        // Ajax request checking state change and replace appropriate element.
        $.ajax(url, {
            type: 'GET',
            dataType: 'json',
            data: {'state': state},
            ContentType: 'application/json',
            success: function (data) {
                // Only replace and activate effects when the state has actually changed.
                if (data.state_changed) {
                    job.find('.state:first').html(data.html);
                    job.transition('pulse');
                    console.log(data.state_changed, state)

                }
                if (data.is_running) {
                    //$("#log").html("");
                    link.attr("href", '/job/view/{0}/#log'.format(uid));
                    loader.html('<div id="log" class="ui log message"><span class="ui active small inline loader"></span><span>Running</span></div>');

                } else {
                    loader.html("");
                    link.attr("href", '/job/view/{0}/'.format(uid))
                }

                imag.replaceWith(data.img_tmpl);
                stdout.text(data.stdout);
                stderr.text(data.stderr);

                // Redirect to the filelist once the job is done.
                if (data.redir && $("#view").length) {
                   window.location.replace(data.redir + "#flist");
                   window.location.reload()
                }

            },
            error: function (xhr, status, text) {
                error_message($(this), xhr, status, text)
            }
        })
    });

}


function preview_template(project_uid) {
    let template = $('#template').val();
    let json_text = $('#json').val();
    $.ajax('/preview/template/',
        {
            type: 'POST',
            dataType: 'json',
            ContentType: 'application/json',
            data: {
                'template': template,
                'json_text': json_text,
                'project_uid': project_uid
            },
            success: function (data) {

                if (data.status === 'success') {
                    $('#preview div').html('<h4 class="ui center aligned header">\n' +
                        '\n' +
                        '    <div class="muted">Press ESC to close window</div>\n' +
                        '</h4>\n' +
                        '\n' +
                        '<pre><code class=" language-bash line-numbers ">' + data.script + '</code></pre>\n');
                    $('#preview').modal('show');
                    Prism.highlightAll();
                } else {
                    popover_message($("#template"), data.msg, data.status);
                }

            },
            error: function (xhr, status, text) {
                error_message($(this), xhr, status, text)
            }

        });
}



function add_vars() {
    let json_text = $('#json').val();
    let template = $('#template').val();

    $.ajax('/add/vars/', {
            type: 'POST',
            dataType: 'json',
            data: {
                'json_text': json_text,
                'template': template,
            },

            success: function (data) {

                $('#template').val(data.code);
                $('#template_field').html(data.html);

            },
            error: function () {
            }
        }
    )
}

function remove_trigger() {
    // Makes site messages dissapear.
    $('.remove').delay(2000).slideUp(800, function () {
        $(this).remove();
    });
}


function set_source_dir() {
    let current_source = $('#current_source');
    if (!current_source.val().length) {
        window.location.href = '/root/list/';
        return
    }
    window.location.href = '/file/list/' + current_source.val();
}


function copy_object(uid, project_uid, clipboard) {
    let container = $('#item-' + uid);

    $.ajax('/copy/object/',
        {
            type: 'POST',
            dataType: 'json',
            data: {
                'project_uid': project_uid,
                'uid': uid,
                'clipboard': clipboard
            },

            success: function (data) {
                if (data.status === 'success') {
                    // Get the item container
                    container.transition({
                        animation: 'pulse', onComplete: function () {
                            container.addClass('copied')
                        }
                    });
                    return
                }
                popover_message(container, data.msg, data.status, 2000)

            },
            error: function (xhr, status, text) {
                error_message($(this), xhr, status, text)
            }


        }
    )
}

function copy_file(path, rel_path) {
    let elem = $('.copy_msg[data-rel="' + rel_path + '"]');
    //alert(elem.html());
    $.ajax('/file/copy/', {
            type: 'POST',
            dataType: 'json',
            data: {
                'path': path
            },

            success: function (data) {
                if (data.status === 'success') {

                    popover_message(elem, data.msg, data.status, 500);
                }
                popover_message(elem, data.msg, data.status, 2000)

            },
            error: function (xhr, status, text) {
                error_message($(this), xhr, status, text)
            }
        }
    )
}


function toggle_delete(elem, otype) {


    var uid = elem.data("value");

    // Fetch the element to decrement/increment  once toggle is complete.
    var count = $("#{0}_count".format(otype));

    var url = '/toggle/delete/';

    var data = { 'uid': uid,  'type': otype };

    $.ajax(url, {
            type: 'POST',
            dataType: 'json',
            data: data,
            success: function (data) {
                if (data.status === 'success') {
                    // Give message
                    elem.transition('zoom');
                    // Update counts on tabular menu
                    count.html(data.counts);
                    return
                }
                popover_message(elem.before(), data.msg, data.status, 2000)

            },
            error: function (xhr, status, text) {
                error_message(elem.before(), xhr, status, text)
            }
        }
    )


}

function change_access(access, user_id, project_uid, elem) {
    let container = $('.container-' + user_id);
    //alert(container.html());
    $.ajax('/manage/access/', {
            type: 'POST',
            dataType: 'json',
            data: {
                'user_id': user_id,
                'access': access,
                'project_uid': project_uid
            },

            success: function (data) {
                if (data.status === 'success') {

                    container.transition('bounce').transition({
                        animation: 'pulse', onComplete: function () {
                            if (data.no_access) {
                                container.removeClass('inputcolor')
                            } else {
                                container.addClass('inputcolor')
                            }
                        }
                    });
                    return
                }
                popover_message(container, data.msg, data.status, 2000)

            },
            error: function (xhr, status, text) {
                error_message(elem, xhr, status, text)
            }
        }
    )

}

$(document).ready(function () {

    $('.access-dropdown').dropdown();
    $(this).on('click', '.access-dropdown .menu .item', function () {
        let user = $(this).parent().data('user');
        let project = $(this).parent().data('project');
        let value = $(this).data("value");
        change_access(value, user, project)

    });

    $('.ui.dropdown').dropdown();

    $('select').dropdown();

    //$('#json_add').dropdown();

    //$('#code_add').dropdown();

    //remove_trigger();

    // Check and update 'Running' and 'Spooled' jobs every 5 seconds.
    setInterval(check_jobs, 5000);

    $(".copy-data").click(function (event) {

        var elem = $(this);
        var data_uid = elem.attr('data-uid');
        var copy_url = elem.attr('copy-url');

        $.ajax(copy_url, {
            type: 'GET',
            dataType: 'json',
            ContentType: 'application/json',
            data: {data_uid: data_uid},
            success: function (data) {
                popover_message($("#copy-message-" + data_uid), data.msg, data.status);
            },
            error: function () {
            }
        });
    });


    $('#recipe-search').keyup(function () {

        var query = $(this).val();

        $.ajax("/search/", {
            type: 'GET',
            dataType: 'html',

            data: {'q': query},

            success: function (data) {
                $('#search-results').html(data);
            },
            error: function () {
            }
        });

    });

    $('#json_add_menu .item').popup({
        on: 'hover'
    });

    $('.listing').popup({
        on: 'hover',
    });


    $('.cmd-value').popup({
        on: 'hover'
    });

    $('#code_add_menu .item').popup({
        on: 'hover'
    });

    $(this).on('click', '#add_vars', function () {
        add_vars()
    });


    $(this).on('click', '#template_preview', function () {
        event.preventDefault();
        //let uid = $(this).data('value');
        let project_uid = $(this).data('value');
        preview_template(project_uid)

    });


    $(this).on('keyup', '#current_source', function (event) {
        // Submit when pressing enter
        if (event.keyCode === 13) {
            set_source_dir()
        }
    });

    $(".moderate-user").click(function (event) {
        event.preventDefault();
        var elem = $(this);

        $('#modpanel').remove();

        // Could be a user or post uid
        var data_uid = elem.attr('data-value');

        var container = $("#moderate-insert-" + data_uid);
        var mod_url = '/accounts/moderate/' + data_uid + '/';

        var page = $('<div id="modpanel"></div>').load(mod_url);
        container.after(page)
    });

    $(this).on('click', '.job .delete', function (event) {
        event.preventDefault();
        var elem = $(this).closest('.job');
        toggle_delete(elem, 'job')
    });

    $(this).on('click', '.data .delete', function (event) {
        event.preventDefault();
        var elem = $(this).closest('.data');
        toggle_delete(elem, 'data')
    });

    $(this).on('click', '#set_source', function () {
        set_source_dir()
    });

    $('.checkbox').checkbox();

    $(this).on('click', '.copy-object', function () {
        let uid = $(this).data('value');
        let project_uid = $(this).data('project');
        let clipboard = $(this).data('clipboard');
        copy_object(uid, project_uid, clipboard);

    });

    $('pre').addClass('language-bash');
    $('code').addClass('language-bash').css('padding', '0');
    Prism.highlightAll();

    $(this).on('click', '.copy_file', function () {
        let path = $(this).data('path');
        let rel_path = $(this).data('rel');
        //alert(rel_path);
        copy_file(path, rel_path)
    });


});
