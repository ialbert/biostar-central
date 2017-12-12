Setup the computer running the biostar engine:

    ansible-playbook setup_computer.yml

Install and link the biostar-engine software:

    ansible-playbook -i hosts setup_engine.yml


Deploy certificates

    sudo certbot --nginx

You should test your configuration at:

* https://www.ssllabs.com/ssltest/analyze.html?d=bioinformatics.recipes
* https://www.ssllabs.com/ssltest/analyze.html?d=data.bioinformatics.recipes
* https://www.ssllabs.com/ssltest/analyze.html?d=www.bioinformatics.recipes
