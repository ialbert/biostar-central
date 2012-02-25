
// comment add
function show_add_comment(parent, post_id){
    csrf_html = parent.find("input[name='csrfmiddlewaretoken']").parent().html()
    $(".comment-form").remove();
    parent.append('\
    <form action="/new/comment/' + post_id + '/" method="post" class="comment-form" id="comment-form">            \
        ' + csrf_html + ' \
        <div class="controls"><textarea class="input-xlarge span8" id="textarea" name="content" rows="3"></textarea></div> \
        <div><a class="btn btn-success" href=\'javascript:document.forms["comment-form"].submit()\'><i class="icon-comment"></i> Add comment</a>          \
        (<span> markdown ok</span>) <a class="btn btn-warning" onclick="javascript:obj=$(\'.comment-form\').remove();"><i class="icon-remove"></i> Cancel</a>   </div>       \
    </form>            \
    '
    )
}

// moderation handler
function mod_link_clicked(link){
    action = link.attr("action")
    table  = link.parents('table') // Find the table holding the entire post to be moderated
    postid = table.find('input[type="hidden"]').val() // Now find the post id from its votebox
    $.post('/moderate/post/' + postid +'/' + action + '/',    
        function(data){
            popover(link.parent(), data.msg, data.status);
            
            function set(value){
                link.text(value); link.attr("action", value)   
            }
            
            if (data.status == 'success') {
                if (action=='delete')  { table.addClass('deleted'); set('undelete') }
                if (action=='undelete'){ table.removeClass('deleted'); set('delete') }
                if (action=='close')   { table.addClass('closed'); set('reopen')  }
                if (action=='reopen')  { table.removeClass('closed'); set('close')  }
            }
        }, 'json');
    
}

//user moderation
function usermod_link_clicked(link){
    info   = $('.user-info')
    action = link.attr("action")
    userid = $('#userid').text() // find the userid
    $.post('/moderate/user/' + userid +'/' + action + '/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            if (data.status == 'success') {
                if (action=='suspend')  { info.addClass('suspended'); }
                if (action=='reinstate')  { info.removeClass('suspended');}
                link.delay(1000).fadeOut(1000, function(){
                    link.remove() 
                });
            }
        }, 'json');
}

//comment deletion
function comment_delete(link){
    pid  = link.attr('target') // post id encoded in link
    par  = $('#'+pid)         // parent div
    body = par.children('div[name="content"]')
    $.post('/post/destroy/' + pid +'/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            if (data.msg == 'destroyed') {
                par.hide('fast')
            } else {
                body.html('[ content removed ]')
            }
        }, 'json');
}
