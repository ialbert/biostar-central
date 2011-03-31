tinyMCE.init({
    theme : "advanced",
    mode : "textareas",
	
	width: "700",
	height: "500",

    //theme_advanced_buttons1 : "bold,italic,underline,undo,redo,link,unlink,image,forecolor,styleselect,removeformat,cleanup,code",
	theme_advanced_buttons1 : "removeformat,styleselect,bold,italic,image,separator,bullist,numlist,separator,link,unlink,separator,code",
    theme_advanced_buttons2 : "",
    theme_advanced_buttons3 : "",
    theme_advanced_toolbar_location : "top",
    theme_advanced_toolbar_align : "center",
    theme_advanced_styles : "Code=codeStyle",
	
	content_css : "/static/biostar.css",
	force_br_newlines : true,
	force_p_newlines : false,
	forced_root_block : '' // Needed for 3.x

	
});