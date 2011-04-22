// Depends on vote.js for popover() function

$(document).ready(function(){
    $('.mod-link').click(function(){
        mod_link_clicked($(this))
    });
});

function mod_link_clicked(link){
    action = link.text()
    table = link.parents('table') // Find the table holding the entire post to be moderated
    post = table.find('input[type="hidden"]').val() // Now find the post id from its votebox
    if(confirm('Are you sure you want to ' + action + ' the post?')){
        $.post('/moderate/' , {post:post, action:action},
        function(data){
            popover(link.parent(), data.msg, data.status)
        }, 'json');
    }
}
