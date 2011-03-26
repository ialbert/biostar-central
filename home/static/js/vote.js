$(document).ready(function(){
    $('.vote').each(function(){
        elem = $(this)
        
        post_id = elem.children('input').val()

        
        up_button = elem.children('.vote-up')
        down_button = elem.children('.vote-down')
        up_button.click(function(){
            do_vote(post_id, 0); 
        });
        down_button.click(function(){
            do_vote(post_id, 1); 
        });
    });
});

function do_vote(post, type){
    //alert('Do vote on post ' + post + ' of type ' + type)
    $.post('/vote/' , {post:post, type:type},
    function(data){
        alert('Received response: ' + data)
    });
}