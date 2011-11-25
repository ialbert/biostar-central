// Depends on vote.js for popover() function

$(document).ready(function(){
    $('.mod-link').click(function(){
        mod_link_clicked($(this))
    });
    $('.usermod-link').click(function(){
        usermod_link_clicked($(this))
    });
    $('.comment-delete').click(function(){
        comment_delete($(this))
    });
});

function mod_link_clicked(link){
    action = link.text()
    table  = link.parents('table') // Find the table holding the entire post to be moderated
    postid = table.find('input[type="hidden"]').val() // Now find the post id from its votebox
    if(confirm('Are you sure you want to ' + action + ' this post')){
        $.post('/moderate/post/' + postid +'/' + action + '/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            window.location.reload()
        }, 'json');
    }
}

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
