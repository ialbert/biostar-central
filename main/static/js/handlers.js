
// comment add
function show_add_comment(parent, post_id){
    csrf_html = parent.find("input[name='csrfmiddlewaretoken']").parent().html()
    parent.html('\
    <form action="/new/comment/' + post_id + '/" method="post" id="comment-form">            \
        ' + csrf_html + ' \
        <div class="controls"><textarea class="input-xlarge span8" id="textarea" name="content" rows="3"></textarea></div> \
        <div><a class="btn btn-success" href=\'javascript:document.forms["comment-form"].submit()\'><i class="icon-comment"></i> Add comment</a>          \
        (<span> markdown ok</span>) <a class="btn btn-warning" href=\'javascript:$("#comment-form")\'><i class="icon-remove"></i> Cancel</a>   </div>       \
    </form>            \
    '    
    )
}

// moderation handler
function mod_link_clicked(link){
    action = link.text()
    table  = link.parents('table') // Find the table holding the entire post to be moderated
    postid = table.find('input[type="hidden"]').val() // Now find the post id from its votebox
    $.post('/moderate/post/' + postid +'/' + action + '/',    
        function(data){
            popover(link.parent(), data.msg, data.status);
            if (action=='delete')  { table.addClass('deleted'); link.text('undelete'); }
            if (action=='undelete'){ table.removeClass('deleted'); link.text('delete'); }
            if (action=='close')   { table.addClass('closed'); link.text('reopen'); }
            if (action=='reopen')  { table.removeClass('closed'); link.text('close'); }
        }, 'json');
    
}

//user moderation
function usermod_link_clicked(link){
    action = link.text()
    userid = $('#userid').text() // find the userid
    username = $('#username').text()
    if(confirm('Are you sure you want to ' + action +' user ' + username + ' ?')){
         $.post('/moderate/user/' + userid +'/' + action + '/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            window.location.reload()
        }, 'json');
    }
}

//comment deletion
function comment_delete(link){
   pid  = link.attr('name')
   if(confirm('Are you sure you want to delete the comment ' + pid + ' ?')){
        $.post('/destroy/post/' + pid +'/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            window.location.reload()
        }, 'json');
    }
}