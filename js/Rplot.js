(function($) {
    $(document).ready(function() {
	
	$('#Rplot').scianimator({
	    'images': ['images/Rplot1.png', 'images/Rplot2.png', 'images/Rplot3.png', 'images/Rplot4.png', 'images/Rplot5.png', 'images/Rplot6.png'],
	    'width': 480,
	    'delay': 500,
	    'loopMode': 'loop'
	});
	$('#Rplot').scianimator('play');
    });
})(jQuery);
