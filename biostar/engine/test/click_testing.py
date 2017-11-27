# -*- coding: utf-8 -*-
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import NoAlertPresentException
import unittest, time, re

class UntitledTestCase(unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(30)
        self.base_url = "https://www.katalon.com/"
        self.verificationErrors = []
        self.accept_next_alert = True
    
    def test_untitled_test_case(self):
        driver = self.driver
        driver.get("http://www.lvh.me:8000/")
        driver.find_element_by_link_text("Project List").click()
        driver.find_element_by_link_text("Login").click()
        driver.find_element_by_id("id_password").clear()
        driver.find_element_by_id("id_password").send_keys("1@lvh.me")
        driver.find_element_by_id("id_email").clear()
        driver.find_element_by_id("id_email").send_keys("1@lvh.me")
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_link_text("Project List").click()
        driver.find_element_by_link_text("Giraffe Genome Assembly and Annotation").click()
        driver.find_element_by_link_text("5 Data Files").click()
        driver.find_element_by_link_text("UGA").click()
        driver.find_element_by_link_text("Edit Data").click()
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_xpath("//a[4]/div").click()
        driver.find_element_by_xpath("//div[4]/div").click()
        driver.find_element_by_xpath("//a[2]/div").click()
        driver.find_element_by_link_text("Analysis Cookbook").click()
        driver.find_element_by_xpath("//div[3]/div/div/div[2]").click()
        driver.find_element_by_link_text("All FASTQ").click()
        driver.find_element_by_link_text("Edit Data").click()
        driver.find_element_by_id("id_sticky").click()
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_xpath("//a[4]/div").click()
        driver.find_element_by_xpath("//a[3]/div").click()
        driver.find_element_by_xpath("//div[3]/div/div[2]").click()
        driver.find_element_by_link_text("FastQC report").click()
        driver.find_element_by_link_text("Recipe Code").click()
        driver.find_element_by_link_text("Edit Script").click()
        driver.find_element_by_name("save_or_preview").click()
        driver.find_element_by_xpath("(//button[@name='save_or_preview'])[2]").click()
        driver.find_element_by_link_text("Back").click()
        driver.find_element_by_xpath("//div[5]/div[3]/a[2]/i").click()
        driver.find_element_by_xpath("//a[contains(text(),'Go to\n                                    Recipe')]").click()
        driver.find_element_by_link_text("Copy Recipe").click()
        driver.find_element_by_xpath("//div[4]/div").click()
        driver.find_element_by_xpath("//a[4]/div").click()
        driver.find_element_by_link_text("Run Analysis").click()
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_xpath("//a[3]/div").click()
        driver.find_element_by_link_text("8 Analysis Results").click()
        driver.find_element_by_xpath("//div[2]/div/div/div[2]").click()
        driver.find_element_by_link_text("Main Result").click()
        driver.find_element_by_link_text("All Generated Files").click()
        driver.find_element_by_xpath("//form/div/div").click()
        driver.find_element_by_xpath("//form/div/div/a/i").click()
        driver.find_element_by_link_text("Go to Recipe").click()
        driver.find_element_by_link_text("Go to Results").click()
        driver.find_element_by_xpath("//a[contains(@href, '/job/view/17/')]").click()
        driver.find_element_by_link_text("Edit Job").click()
        driver.find_element_by_id("id_name").click()
        driver.find_element_by_id("id_name").send_keys(Keys.DOWN)
        driver.find_element_by_id("id_name").clear()
        driver.find_element_by_id("id_name").send_keys("FastQC report2")
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_xpath("//a[4]/div").click()
        driver.find_element_by_xpath("//a[3]/div").click()
        driver.find_element_by_link_text("Edit Project").click()
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_xpath("//div[3]/div/div[4]").click()
        driver.find_element_by_xpath("//a[2]/div").click()
        driver.find_element_by_link_text("Analysis Cookbook").click()
        driver.find_element_by_link_text("Manage People").click()
        driver.find_element_by_name("searches").click()
        driver.find_element_by_name("searches").clear()
        driver.find_element_by_name("searches").send_keys("a")
        driver.find_element_by_xpath("//button/i").click()
        driver.find_element_by_name("users").click()
        driver.find_element_by_name("add_or_remove").click()
        driver.find_element_by_xpath("//a[3]/div").click()
        driver.find_element_by_xpath("//a[2]/div").click()
        driver.find_element_by_link_text("1@lvh.me").click()
        driver.find_element_by_link_text("Edit profile").click()
        driver.find_element_by_id("id_first_name").click()
        driver.find_element_by_link_text("Reset password").click()
        driver.find_element_by_xpath("//a/div").click()
        driver.find_element_by_link_text("1@lvh.me").click()
        driver.find_element_by_link_text("Edit profile").click()
        driver.find_element_by_id("id_first_name").click()
        driver.find_element_by_id("id_first_name").clear()
        driver.find_element_by_id("id_first_name").send_keys("Admin Userr")
        driver.find_element_by_xpath("//button[@type='submit']").click()
        driver.find_element_by_link_text("Home").click()
        driver.find_element_by_link_text("Admin").click()
        driver.find_element_by_link_text("Logout").click()
        driver.find_element_by_xpath("//button[@type='submit']").click()
    
    def is_element_present(self, how, what):
        try: self.driver.find_element(by=how, value=what)
        except NoSuchElementException as e: return False
        return True
    
    def is_alert_present(self):
        try: self.driver.switch_to_alert()
        except NoAlertPresentException as e: return False
        return True
    
    def close_alert_and_get_its_text(self):
        try:
            alert = self.driver.switch_to_alert()
            alert_text = alert.text
            if self.accept_next_alert:
                alert.accept()
            else:
                alert.dismiss()
            return alert_text
        finally: self.accept_next_alert = True
    
    def tearDown(self):
        self.driver.quit()
        self.assertEqual([], self.verificationErrors)

if __name__ == "__main__":
    unittest.main()
