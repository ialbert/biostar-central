// Depends on vote.js for popover() function

function comment_delete(link){
   pid  = link.attr('name')
   if(confirm('Are you sure you want to delete the comment ' + pid + ' ?')){
        $.post('/comment/delete/' + pid +'/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            window.location.reload()
        }, 'json');
    }
}
