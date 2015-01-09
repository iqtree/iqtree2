/*
 * modelsblock.cpp
 *
 *  Created on: Jan 9, 2015
 *      Author: minh
 */

#include "modelsblock.h"

ModelsBlock::ModelsBlock()  : NxsBlock()
{
	id = "MODELS";
}

ModelsBlock::~ModelsBlock() {
}

void ModelsBlock::Read(NxsToken &token)
{
	// This should be the semicolon after the block name
	token.GetNextToken();

	if (!token.Equals(";"))
		throw NxsException("Expecting ';' after MODELS block name", token);
	for (;;) {
		token.GetNextToken();
		if (token.Equals("MODEL")) {
			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextToken();
			NxsModel model;

			model.name = token.GetToken();

			token.GetNextToken();
			if (!token.Equals("="))
				throw NxsException("Expecting '=' after model name", token);

			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextContiguousToken(';');
			model.description = token.GetToken();

			token.GetNextToken();
			if (!token.Equals(";"))
				throw NxsException("Expecting ';' to terminate MODEL command", token);

			push_back(model);

		} else if (token.Equals("END") || token.Equals("ENDBLOCK")) {
			// Get the semicolon following END
			token.GetNextToken();

			if (!token.Equals(";"))
				throw NxsException("Expecting ';' to terminate the ENDBLOCK command, but found ", token);
			break;
		}	// if (token.Equals("END") || token.Equals("ENDBLOCK"))
		else {
			SkippingCommand(token.GetToken());
			do {
				token.GetNextToken();
			} while (!token.AtEOF() && !token.Equals(";"));

			if (token.AtEOF())
				throw NxsException("Unexpected end of file encountered", token);
		}	// token not END, ENDBLOCK, COST
	}
}

NxsModel *ModelsBlock::findModel(string name) {
	for (iterator it = begin(); it != end(); it++)
		if (it->name == name) return &(*it);
	return NULL;
}
