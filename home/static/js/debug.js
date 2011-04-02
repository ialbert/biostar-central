// Sets up jQuery to pop-up the full Django error log in case of a
// 500 reponse to an AJAX request

$('html').ajaxError(function(event, xhr, settings, error){
    //elem =
    elem = $('body').after('<div class="popup"></div>')
    $('.popup').css({'position':'absolute', 'background-color':'#FFFFCC', 'top':100, 'left':100, 'width':'80%', 'z-index':1});
    $('.popup').html(xhr.responseText);
    $('.popup style').remove();
    $('.popup').prepend('<div class="close-button">Close</div>');
    $('.close-button').css({'float':'right'});
    $('.close-button').click(function(){
       $('.popup').remove(); 
    });
    $('.popup').dblclick(function(){
       $(this).remove(); 
        
    });
    //alert(xhr.responseText);
    //newwindow=window.open(url,'name','height=400,width=200');
});

// Also controls the SQL queries item in the footer in debug mode

$(document).ready(function(){
    $('#queries').hide()
    $('#toggle-queries').click(function(){
       $('#queries').toggle() 
    });
});