set(TEST_FILES
	test1.sh
	test2.sh
	test3.sh
	)

foreach(TEST_FILE ${TEST_FILES})
	add_test(NAME ${TEST_FILE} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMAND sh ./tests/${TEST_FILE})
endforeach()
